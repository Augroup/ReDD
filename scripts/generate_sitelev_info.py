from pathlib import Path
import multiprocessing as mp
import time
import argparse
import pickle
import shutil
from patch_mp import patch_mp_connection_bpo_17560
patch_mp_connection_bpo_17560()
def get_file_marker(filename,threads):
    file_size = Path(filename).stat().st_size
    single_thread_file_size, extra = divmod(file_size, threads)
    if extra:
        single_thread_file_size += 1
    list_of_start_pos = []
    with open(filename) as f:
        for t in range(threads):
            if t == 0:
                start_pos = 0
                list_of_start_pos.append(start_pos)
            else:
                f.seek(t*single_thread_file_size)
                start_pos = len(f.readline().encode('utf8')) + t*single_thread_file_size
                list_of_start_pos.append(start_pos)
    return list_of_start_pos,file_size
import hashlib
def chr_name_hash(chr_name,num_bins):
    first_byte = hashlib.md5(chr_name.encode('utf8')).digest()[0]
    second_byte = hashlib.md5(chr_name.encode('utf8')).digest()[1]
    return (first_byte*10+second_byte)%num_bins
def dump_info(pos_coverage,predict_scores,coverage,strands,temp_folder,worker_id,hashed_chr_name_set,num_bins):
    for chr_name in pos_coverage:
        hashed_chr_name = chr_name_hash(chr_name,num_bins)
        hashed_chr_name_set.add(hashed_chr_name)
        with open(f'{temp_folder}/{hashed_chr_name}_{worker_id}','ab') as f:
            pickle.dump(({chr_name:pos_coverage[chr_name]},{chr_name:predict_scores[chr_name]},{chr_name:coverage[chr_name]},{chr_name:strands[chr_name]}),f)
    return hashed_chr_name_set
    # with open(f'{temp_folder}/{worker_id}','ab') as f:
    #     pickle.dump((pos_coverage,predict_scores,coverage,strands),f)
# def agg_info(temp_folder,worker_id,num_bins):
#     all_pos_coverage={}
#     all_predict_scores={}
#     all_coverage = {}
#     all_strands ={}
#     with open(f'{temp_folder}/{worker_id}','rb') as f:
#         while 1:
#             try:
#                 (pos_coverage,predict_scores,coverage,strands) = pickle.load(f)
#                 for chr_name in pos_coverage:
#                     if chr_name not in all_pos_coverage:
#                         all_pos_coverage[chr_name] = {}
#                         all_coverage[chr_name] = {}
#                         all_predict_scores[chr_name]={}
#                         all_strands[chr_name]={}
#                     for pos in pos_coverage[chr_name]:
#                         if pos not in all_pos_coverage[chr_name]:
#                             all_pos_coverage[chr_name][pos] = 0
#                             all_predict_scores[chr_name][pos]=0
#                             all_coverage[chr_name][pos] = 0
#                         all_pos_coverage[chr_name][pos] += pos_coverage[chr_name][pos]
#                         all_coverage[chr_name][pos] += coverage[chr_name][pos]
#                         all_predict_scores[chr_name][pos] += predict_scores[chr_name][pos]
#                         all_strands[chr_name][pos] = strands[chr_name][pos]
#             except EOFError:
#                 break
#     hashed_chr_name_set = set()
#     for chr_name in all_pos_coverage:
#         hashed_chr_name = chr_name_hash(chr_name,num_bins)
#         hashed_chr_name_set.add(hashed_chr_name)
#         with open(f'{temp_folder}/{hashed_chr_name}_{worker_id}','ab') as f:
#             pickle.dump(({chr_name:all_pos_coverage[chr_name]},{chr_name:all_predict_scores[chr_name]},{chr_name:all_coverage[chr_name]},{chr_name:all_strands[chr_name]}),f)
#     return hashed_chr_name_set
def process_coverage(start_pos,end_pos,cutoff,filename,temp_folder,worker_id,num_bins):
    pos_coverage={}
    predict_scores={}
    coverage = {}
    strands ={}
    hashed_chr_name_set = set()
    buffer_size = 2e7
    with open(filename) as f:
        f.seek(start_pos)
        processed_bytes = 0
        processed_lines = 0
        for line in f:
            if len(line.split("\t"))<5:
                continue
            
            score=float(line.split("\t")[-1])
            strand=line.split("\t")[-2]
            processed_bytes += len(line.encode('utf8'))
            if processed_bytes > end_pos - start_pos:
                break
            processed_lines += 1
            if score>cutoff:
                predict_label=1
            else:
                predict_label=0
            
            chr_name= line.split("\t")[2]
            pos = int(line.split("\t")[3])
            if chr_name not in coverage:
                pos_coverage[chr_name] = {}
                coverage[chr_name] = {}
                predict_scores[chr_name] = {}
                strands[chr_name]={}
            if pos not in coverage[chr_name]:
                pos_coverage[chr_name][pos]=predict_label
                predict_scores[chr_name][pos]=score
                coverage[chr_name][pos]=1
                strands[chr_name][pos] = strand
            else:
                pos_coverage[chr_name][pos]+=predict_label
                predict_scores[chr_name][pos]+=score
                coverage[chr_name][pos]+=1
            if processed_lines >= buffer_size:
                hashed_chr_name_set = dump_info(pos_coverage,predict_scores,coverage,strands,temp_folder,worker_id,hashed_chr_name_set,num_bins)
                processed_lines = 0
                pos_coverage={}
                predict_scores={}
                coverage = {}
                strands ={}
    if processed_lines != 0:
        hashed_chr_name_set = dump_info(pos_coverage,predict_scores,coverage,strands,temp_folder,worker_id,hashed_chr_name_set,num_bins)
        processed_lines = 0
        pos_coverage={}
        predict_scores={}
        coverage = {}
        strands ={}
    # hashed_chr_name_set = agg_info(temp_folder,worker_id,num_bins)
    return {worker_id:hashed_chr_name_set}
def agg_hashed(reverse_all_hashed_chr_name_dict,worker_hashed_chr_name_set,cov_cutoff,ratio_cutoff,mod_cov_cutoff):
    for hashed_chr_name in worker_hashed_chr_name_set:
        all_pos_coverage = {}
        all_coverage = {}
        all_predict_scores={}
        all_strands= {}
        for worker_id in reverse_all_hashed_chr_name_dict[hashed_chr_name]:
            with open(f'{temp_folder}/{hashed_chr_name}_{worker_id}','rb') as f:
                     while 1:
                        try:
                            (pos_coverage,predict_scores,coverage,strands) = pickle.load(f)
                            for chr_name in pos_coverage:
                                if chr_name not in all_pos_coverage:
                                    all_pos_coverage[chr_name] = {}
                                    all_coverage[chr_name] = {}
                                    all_predict_scores[chr_name]={}
                                    all_strands[chr_name]={}
                                for pos in pos_coverage[chr_name]:
                                    if pos not in all_pos_coverage[chr_name]:
                                        all_pos_coverage[chr_name][pos] = 0
                                        all_predict_scores[chr_name][pos]=0
                                        all_coverage[chr_name][pos] = 0
                                    all_pos_coverage[chr_name][pos] += pos_coverage[chr_name][pos]
                                    all_coverage[chr_name][pos] += coverage[chr_name][pos]
                                    all_predict_scores[chr_name][pos] += predict_scores[chr_name][pos]
                                    all_strands[chr_name][pos] = strands[chr_name][pos]
                        except EOFError:
                            break
        output_dict = {}
        for chr_name in all_coverage:
            for pos in sorted(all_coverage[chr_name].keys()):
                ratio = float(all_pos_coverage[chr_name][pos]/all_coverage[chr_name][pos])
                avg_predict_score=float(all_predict_scores[chr_name][pos]/all_coverage[chr_name][pos])
                if cov_cutoff is not None:
                    if all_coverage[chr_name][pos] < cov_cutoff:
                        continue
                if ratio_cutoff is not None:
                    if ratio < ratio_cutoff:
                        continue
                if mod_cov_cutoff is not None:
                    if all_pos_coverage[chr_name][pos] < mod_cov_cutoff:
                        continue
                if chr_name not in output_dict:
                    output_dict[chr_name] = {}
                if pos not in output_dict[chr_name]:
                    output_dict[chr_name][pos] = {}
                output_dict[chr_name][pos] = {'pos_coverage':all_pos_coverage[chr_name][pos],
                'coverage':all_coverage[chr_name][pos],'strands':all_strands[chr_name][pos],
                'ratio':ratio,'avg_predict_score':avg_predict_score}
        with open(f'{temp_folder}/{hashed_chr_name}','wb') as f:
            pickle.dump(output_dict,f)


def aggregate_and_generate_output(threads,temp_folder,outputname,reverse_all_hashed_chr_name_dict):
        for hashed_chr_name in reverse_all_hashed_chr_name_dict:
            # all_pos_coverage = {}
            # all_coverage = {}
            # all_predict_scores={}
            # all_strands= {}
            # for worker_id in reverse_all_hashed_chr_name_dict[hashed_chr_name]:
            #     with open(f'{temp_folder}/{hashed_chr_name}_{worker_id}','rb') as f:
            #          while 1:
            #             try:
            #                 (pos_coverage,predict_scores,coverage,strands) = pickle.load(f)
            #                 for chr_name in pos_coverage:
            #                     if chr_name not in all_pos_coverage:
            #                         all_pos_coverage[chr_name] = {}
            #                         all_coverage[chr_name] = {}
            #                         all_predict_scores[chr_name]={}
            #                         all_strands[chr_name]={}
            #                     for pos in pos_coverage[chr_name]:
            #                         if pos not in all_pos_coverage[chr_name]:
            #                             all_pos_coverage[chr_name][pos] = 0
            #                             all_predict_scores[chr_name][pos]=0
            #                             all_coverage[chr_name][pos] = 0
            #                         all_pos_coverage[chr_name][pos] += pos_coverage[chr_name][pos]
            #                         all_coverage[chr_name][pos] += coverage[chr_name][pos]
            #                         all_predict_scores[chr_name][pos] += predict_scores[chr_name][pos]
            #                         all_strands[chr_name][pos] = strands[chr_name][pos]
            #             except EOFError:
            #                 break
            with open(f'{temp_folder}/{hashed_chr_name}','rb') as f:
                output_dict = pickle.load(f)
                with open(outputname,'w') as output:
                    for chr_name in output_dict:
                        for pos in sorted(output_dict[chr_name].keys()):
                            pos_coverage_val = output_dict[chr_name][pos]['pos_coverage']
                            coverage_val = output_dict[chr_name][pos]['coverage']
                            strands_val = output_dict[chr_name][pos]['strands']
                            ratio =  output_dict[chr_name][pos]['ratio']
                            avg_predict_score =  output_dict[chr_name][pos]['avg_predict_score']
                            output.write(f'{chr_name}\t{pos}\t{pos_coverage_val}\t{coverage_val}\t{ratio:.3f}\t{avg_predict_score:.3f}\t{strands_val}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default=None,required=True,help="input")
    parser.add_argument('--output', type=str, default=None,required=True,help="output")
    parser.add_argument('--threads', type=int, default=1,required=False,help="threads")
    parser.add_argument('--cov_cutoff', type=float, default=None,required=False,help="read coverage cutoff")
    parser.add_argument('--ratio_cutoff', type=float, default=None,required=False,help="site-level A2G ratio cutoff")
    parser.add_argument('--mod_cov_cutoff', type=float, default=None,required=False,help="modified read coverage cutoff")
    parser.add_argument('--temp_dir', type=float, default=None,required=False,help="temp dir")
    args = parser.parse_args()
    for k, v in vars(args).items():
        print(k, ':', v)
    
    filename=args.input
    outputname=args.output
    cutoff = 0.5
    threads = args.threads
    cov_cutoff=args.cov_cutoff
    ratio_cutoff=args.ratio_cutoff
    mod_cov_cutoff=args.mod_cov_cutoff
    temp_folder = args.temp_dir
    if temp_folder is None:
        if '/' in outputname:
            temp_folder = '/'.join(outputname.split('/')[:-1])+'/temp/'+outputname.split('/')[-1].split('.')[0]
        else:
            temp_folder = './temp/'+outputname.split('.')[0]
    if Path(temp_folder).exists():
        print('temp folder is not empty! Will completely remove all the files inside!',flush=True)
        shutil.rmtree(temp_folder, ignore_errors=True)
    Path(temp_folder).mkdir(exist_ok=True,parents=True)
    
    start_time = time.time()
    list_of_start_pos,file_size = get_file_marker(filename,threads)
    num_bins = max(int(file_size // (1024**3)),1)
    print(num_bins,flush=True)
    if threads == 1:
        list_of_end_pos = [file_size]
    else:
        list_of_end_pos = [pos for pos in list_of_start_pos[1:]]
        list_of_end_pos.append(file_size)
    pool = mp.Pool(threads)
    futures = []
    for start_pos,end_pos,worker_id in zip(list_of_start_pos,list_of_end_pos,range(threads)):
        futures.append(pool.apply_async(process_coverage,(start_pos,end_pos,cutoff,filename,temp_folder,worker_id,num_bins)))
    all_hashed_chr_name_dict = {}
    for future in futures:
        hashed_chr_name_dict = future.get()
        all_hashed_chr_name_dict.update(hashed_chr_name_dict)
    pool.close()
    pool.join()
    reverse_all_hashed_chr_name_dict = {}
    for worker_id in all_hashed_chr_name_dict:
        for hashed_chr_name in all_hashed_chr_name_dict[worker_id]:
            if hashed_chr_name not in reverse_all_hashed_chr_name_dict:
                reverse_all_hashed_chr_name_dict[hashed_chr_name] = set()
            reverse_all_hashed_chr_name_dict[hashed_chr_name].add(worker_id)
    all_hashed_chr_name_list = list(set(reverse_all_hashed_chr_name_dict.keys()))
    chunksize,extra = divmod(len(all_hashed_chr_name_list), threads)
    if extra:
        chunksize += 1
    pool = mp.Pool(threads)
    futures = []
    for i in range(threads):
        if (i+1)*chunksize > len(all_hashed_chr_name_list):
            worker_hashed_chr_name_set = all_hashed_chr_name_list[i*chunksize:]
        else:
            worker_hashed_chr_name_set = all_hashed_chr_name_list[i*chunksize:(i+1)*chunksize]
        futures.append(pool.apply_async(agg_hashed,(reverse_all_hashed_chr_name_dict,worker_hashed_chr_name_set,cov_cutoff,ratio_cutoff,mod_cov_cutoff)))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
    aggregate_and_generate_output(threads,temp_folder,outputname,reverse_all_hashed_chr_name_dict)
    running_time = time.time() - start_time
    print(f"Done in {running_time} s")
    # shutil.rmtree(temp_folder, ignore_errors=True)



# output.close()
#datatype="NA12878_dRNA"
#folder="/fs/ess/scratch/PCON0009/duolin/modification/Eventfeature_5feature_win9/merge_8alldata2_del_H1/"
#input=${folder}${datatype}"_noncandidate.txt" 
#output=${folder}${datatype}".sitelev.bed"
#python3 generate_sitelev_info.py --input $input --output $output --threads