from pathlib import Path
import multiprocessing as mp
import time
import argparse

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

def process_coverage(start_pos,end_pos,cutoff,filename):
    pos_coverage={}
    predict_scores={}
    coverage = {}
    with open(filename) as f:
        f.seek(start_pos)
        processed_bytes = 0
        for line in f:
            if len(line.split("\t"))<5:
                continue
            
            score=float(line.split("\t")[-1])
            processed_bytes += len(line.encode('utf8'))
            if processed_bytes > end_pos - start_pos:
                break
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
            if pos not in coverage[chr_name]:
                pos_coverage[chr_name][pos]=predict_label
                predict_scores[chr_name][pos]=score
                coverage[chr_name][pos]=1
            else:
                pos_coverage[chr_name][pos]+=predict_label
                predict_scores[chr_name][pos]+=score
                coverage[chr_name][pos]+=1
    # buf_size = 4 * 1024 * 1024 * 1024
    # if buf_size > end_pos - start_pos:
    #     buf_size = end_pos - start_pos
    # with open(filename) as f:
    #     f.seek(start_pos)
    #     processed_bytes = 0
    #     is_finished = False
    #     while not is_finished:
    #         lines = f.readlines(buf_size)
    #         if len(lines) == 0:
    #             break
    #         for line in lines:
    #             processed_bytes += len(line.encode('utf8'))
    #             if processed_bytes > end_pos - start_pos:
    #                 is_finished = True
    #                 break
    #             if float(line.split("\t")[-1])>cutoff:
    #                 predict_label=1
    #             else:
    #                 predict_label=0
                
    #             chr_name= line.split("\t")[2]
    #             pos = int(line.split("\t")[3])
    #             if chr_name not in coverage:
    #                 pos_coverage[chr_name] = {}
    #                 coverage[chr_name] = {}
    #             if pos not in coverage[chr_name]:
    #                 pos_coverage[chr_name][pos]=predict_label
    #                 coverage[chr_name][pos]=1
    #             else:
    #                 pos_coverage[chr_name][pos]+=predict_label
    #                 coverage[chr_name][pos]+=1
    return pos_coverage,coverage,predict_scores
def aggregate_and_generate_output(list_of_pos_coverage,list_of_coverage,list_of_predict_scores,outputname,cov_cutoff,ratio_cutoff,mod_cov_cutoff):
    all_pos_coverage = list_of_pos_coverage[0]
    all_coverage = list_of_coverage[0]
    all_predict_scores=list_of_predict_scores[0]
    if len(list_of_pos_coverage) != 1:
        for pos_coverage,coverage,predict_scores in zip(list_of_pos_coverage[1:],list_of_coverage[1:],list_of_predict_scores[1:]):
            for chr_name in pos_coverage:
                if chr_name not in all_pos_coverage:
                    all_pos_coverage[chr_name] = {}
                    all_coverage[chr_name] = {}
                    all_predict_scores[chr_name]={}
                for pos in pos_coverage[chr_name]:
                    if pos not in all_pos_coverage[chr_name]:
                        all_pos_coverage[chr_name][pos] = 0
                        all_predict_scores[chr_name][pos]=0
                        all_coverage[chr_name][pos] = 0
                    all_pos_coverage[chr_name][pos] += pos_coverage[chr_name][pos]
                    all_coverage[chr_name][pos] += coverage[chr_name][pos]
                    all_predict_scores[chr_name][pos] += predict_scores[chr_name][pos]
    with open(outputname,'w') as output:
        for chr_name in sorted(all_coverage.keys()):
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
                
                output.write(f'{chr_name}\t{pos}\t{all_pos_coverage[chr_name][pos]}\t{all_coverage[chr_name][pos]}\t{ratio:.3f}\t{avg_predict_score:.3f}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default=None,required=True,help="input")
    parser.add_argument('--output', type=str, default=None,required=True,help="output")
    parser.add_argument('--threads', type=int, default=1,required=False,help="threads")
    parser.add_argument('--cov_cutoff', type=float, default=None,required=False,help="read coverage cutoff")
    parser.add_argument('--ratio_cutoff', type=float, default=None,required=False,help="site-level A2G ratio cutoff")
    parser.add_argument('--mod_cov_cutoff', type=float, default=None,required=False,help="modified read coverage cutoff")
    
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
    start_time = time.time()
    list_of_start_pos,file_size = get_file_marker(filename,threads)
    if threads == 1:
        list_of_end_pos = [file_size]
    else:
        list_of_end_pos = [pos for pos in list_of_start_pos[1:]]
        list_of_end_pos.append(file_size)
    pool = mp.Pool(threads)
    futures = []
    for start_pos,end_pos in zip(list_of_start_pos,list_of_end_pos):
        futures.append(pool.apply_async(process_coverage,(start_pos,end_pos,cutoff,filename)))
    list_of_pos_coverage,list_of_coverage,list_of_predict_scores = [],[],[]
    for future in futures:
        pos_coverage,coverage,predict_scores = future.get()
        list_of_pos_coverage.append(pos_coverage)
        list_of_coverage.append(coverage)
        list_of_predict_scores.append(predict_scores)
    pool.close()
    pool.join()
    
    aggregate_and_generate_output(list_of_pos_coverage,list_of_coverage,list_of_predict_scores,outputname,cov_cutoff,ratio_cutoff,mod_cov_cutoff)
    running_time = time.time() - start_time
    print(f"Done in {running_time} s")



# output.close()
#datatype="NA12878_dRNA"
#folder="/fs/ess/scratch/PCON0009/duolin/modification/Eventfeature_5feature_win9/merge_8alldata2_del_H1/"
#input=${folder}${datatype}"_noncandidate.txt" 
#output=${folder}${datatype}".sitelev.bed"
#python3 generate_sitelev_info.py --input $input --output $output --threads