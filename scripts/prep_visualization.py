import pandas as pd
import pyBigWig
import pysam
import numpy as np
from pathlib import Path
import argparse
import time
import multiprocessing as mp
import os
import glob
def get_color(score):
    score = float(score)
    r_val = int((1-score) * 255)
    b_val = int(score * 255)
    return f"{r_val},0,{b_val}"
def get_cigar(sites_list,score_list,thres):
    cigar_list = []
    reference_start = sites_list[0]
    reference_start_score = score_list[0]
    previous_pos = reference_start
    if reference_start_score >= thres:
        cigar_list.append((0,1))
    else:
        cigar_list.append((2,1))
    if len(sites_list) > 1:
        for pos,score in zip(sites_list[1:],score_list[1:]):
            distance = pos - previous_pos + 1 - 2
            if distance > 0:
                cigar_list.append((3,distance))
                if score >= thres:
                    cigar_list.append((0,1))
                else:
                    cigar_list.append((2,1))
            elif distance == 0:
                if score >= thres:
                    cigar_list.append((0,1))
                else:
                    cigar_list.append((2,1))
            previous_pos = pos
    return cigar_list
def get_line_marker(mole_path,threads):
    with open(mole_path, 'r') as f:
        line_offset = []
        offset = 0
        old_read_id = ''
        is_first_read = True
        for line in f:
            read_id = line.split('\t')[1]
            if read_id != old_read_id and not is_first_read:
                line_offset.append(offset)
                old_read_id = read_id
            is_first_read = False
    #         line_offset.append(offset)
            offset += len(line)
        num_reads = len(line_offset)
        print(f'Number of reads is {num_reads},',flush=True)
        chunksize, extra = divmod(num_reads, threads)
        if extra:
            chunksize += 1
        line_marker = []
        for i in range(threads):
            line_marker.append((line_offset[i*chunksize],chunksize))
    return line_marker
def mole_info_to_bam_worker(mole_path,cDNA_mole_input,reference_id_dict,reference_path,mole_output_dir,i,thres,marker):
    start_file_pos,num_assigned_reads = marker
    read_isoform_mapping = {}
    if cDNA_mole_input is not None:
        with open(cDNA_mole_input) as f:
            for line in f:
                fields = line.split('\t')
                read_id = fields[1]
                isoform = fields[2]
                read_isoform_mapping[read_id] = isoform
    reference_fasta = pysam.FastaFile(reference_path)
    num_processsed_reads = 0
    old_read_id = ''
    old_read_sites = []
    old_fields = []
    temp_folder = f'{mole_output_dir}/temp/'
    Path(temp_folder).mkdir(exist_ok=True,parents=True)
    output_handle = pysam.AlignmentFile(f'{temp_folder}/{i}.bam', "wb",reference_names=reference_fasta.references,reference_lengths=reference_fasta.lengths,threads=1)
    with open(mole_path,'r') as f:
        f.seek(start_file_pos)
        for line in f:
            fields = line.split('\t')
            [_,read_id,chr_name,start_pos,strand,score] = fields
#             if float(score) < thres:
#                 continue
            if chr_name not in reference_id_dict:
                continue
            start_pos = int(start_pos)
            if read_id != old_read_id and old_read_id != '':
                if num_processsed_reads > num_assigned_reads:
                    print(f'Process {i} has handled {num_processsed_reads} reads!',flush=True)
                    break
                num_processsed_reads += 1
                score_list = []
                sites_list = []
                for pos,score in old_read_sites:
                    score_list.append(float(score))
                    sites_list.append(pos)
                seg = pysam.AlignedSegment()
                seg.query_name = str(old_read_id)
                seg.flag = 0 if old_fields[4] == '+' else 16
                seg.cigar = get_cigar(sites_list,score_list,thres)
                seg.reference_id = reference_id_dict[old_fields[2]]
                seg.reference_start = sites_list[0] - 1
                seg.mapping_quality = 60
                score_str = str(round(np.mean(score_list),2))
#                 seg.tags = [('SC',get_color(np.mean(score_list)),'Z'),('ML',array.array('f',score_list),'f')]
                tags = [('SC',get_color(np.mean(score_list)),'Z'),('ML',f'Score:{score_str}','Z')]
                if str(old_read_id) in read_isoform_mapping:
                    isoform = read_isoform_mapping[old_read_id]
                    tags.append(('TT',isoform,'Z'))
                seg.tags = tags
                output_handle.write(seg)
                # new read
                old_read_id = read_id
                old_read_sites = []
            old_read_id = read_id
            old_read_sites.append((start_pos,score))   
            old_fields = fields
        # if reach last line of file
        if len(old_fields) != 0:
            score_list = []
            sites_list = []
            for pos,score in old_read_sites:
                score_list.append(float(score))
                sites_list.append(pos)
            seg = pysam.AlignedSegment()
            seg.query_name = str(old_read_id)
            seg.flag = 0 if old_fields[4] == '+' else 16
            seg.cigar = get_cigar(sites_list,score_list,thres)
            seg.reference_id = reference_id_dict[old_fields[2]]
            seg.reference_start = sites_list[0] - 1
            seg.mapping_quality = 60
            score_str = str(round(np.mean(score_list),2))
#                 seg.tags = [('SC',get_color(np.mean(score_list)),'Z'),('ML',array.array('f',score_list),'f')]
            tags = [('SC',get_color(np.mean(score_list)),'Z'),('ML',f'Score:{score_str}','Z')]
            if str(old_read_id) in read_isoform_mapping:
                isoform = read_isoform_mapping[old_read_id]
                tags.append(('TT',isoform,'Z'))
            seg.tags = tags
            output_handle.write(seg)
            output_handle.close()
def mole_info_to_bam(mole_path,reference_id_dict,reference_path,mole_output_dir,cDNA_mole_input,threads):
    line_marker = get_line_marker(mole_path,threads)
    temp_folder = f'{mole_output_dir}/temp/'
    for thres in range(10):
        thres = thres / 10
        pool = mp.Pool(threads+1)
        futures = []
        all_bam_files = []
        for marker,i in zip(line_marker,range(threads)):
            all_bam_files.append(f'{temp_folder}/{i+1}.bam')
            futures.append(pool.apply_async(mole_info_to_bam_worker,(mole_path,cDNA_mole_input,reference_id_dict,reference_path,mole_output_dir,i+1,thres,marker,)))
        for future in futures:
            future.get()
        sorted_bam = f'{mole_output_dir}/molecule_thres_{int(thres*10)}.sorted.bam'
        bam = f'{mole_output_dir}/molecule_thres_{int(thres*10)}.bam'
        params = ['-@',str(threads),'-f',bam]+all_bam_files
        pysam.merge(*params)
        pysam.sort('-@',str(threads),"-o",sorted_bam,bam)
        pysam.index(sorted_bam)
        for f in glob.glob(os.path.join(temp_folder, "*.bam")):
            os.remove(f)
    os.rmdir(temp_folder)
def pre_compute_site_level_thres(output_dir,header,sitelevel_df,candidate_df):
    sitelevel_df = sitelevel_df.sort_values(by=['rname','pos'])
    if candidate_df is not None:
        indexed_candidate_df = candidate_df.set_index(['rname','pos'])
        indexed_sitelevel_df = sitelevel_df.set_index(['rname','pos'])
        candidate_idx = indexed_sitelevel_df.index.intersection(indexed_candidate_df.index)
        for thres in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
            bw = pyBigWig.open(f"{output_dir}/cov_score_candidate_thres_{int(thres*10)}.bw", "wb")
            bw.addHeader(header,maxZooms=0)
            candidate_sitelevel_df = indexed_sitelevel_df.loc[candidate_idx].reset_index()
            candidate_sitelevel_df = candidate_sitelevel_df[(candidate_sitelevel_df['cov_score'] >= thres)].sort_values(by=['rname','pos'])
            # index = sitelevel_df.index
            chroms = candidate_sitelevel_df['rname'].values
            starts = np.array(candidate_sitelevel_df['start'].values, dtype=np.int64)
            ends = np.array(candidate_sitelevel_df['end'], dtype=np.int64)
            score = np.round(np.array(candidate_sitelevel_df['cov_score'], dtype=np.float),4)

            bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=score.tolist(),validate=True)

            bw.close()
        for thres in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
            bw = pyBigWig.open(f"{output_dir}/cov_score_noncandidate_thres_{int(thres*10)}.bw", "wb")

            bw.addHeader(header,maxZooms=0)
            non_candidate_sitelevel_df = indexed_sitelevel_df.loc[~indexed_sitelevel_df.index.isin(candidate_idx)].reset_index()
            non_candidate_sitelevel_df = non_candidate_sitelevel_df[(non_candidate_sitelevel_df['cov_score'] >= thres)].sort_values(by=['rname','pos'])

            # index = sitelevel_df.index
            chroms = non_candidate_sitelevel_df['rname'].values
            starts = np.array(non_candidate_sitelevel_df['start'].values, dtype=np.int64)
            ends = np.array(non_candidate_sitelevel_df['end'], dtype=np.int64)
            score = np.round(np.array(non_candidate_sitelevel_df['cov_score'], dtype=np.float),4)

            bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=score.tolist(),validate=True)

            bw.close()
    else:  
        for thres in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
            bw = pyBigWig.open(f"{output_dir}/cov_score_thres_{int(thres*10)}.bw", "wb")

            bw.addHeader(header,maxZooms=0)
            index = (sitelevel_df['cov_score'] >= thres)
            chroms = sitelevel_df.loc[index,'rname'].values
            starts = np.array(sitelevel_df.loc[index,'start'].values, dtype=np.int64)
            ends = np.array(sitelevel_df.loc[index,'end'], dtype=np.int64)
            score = np.array(sitelevel_df.loc[index,'cov_score'])

            bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=score.tolist(),validate=True)

            bw.close()
    for thres in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
        not_pass_sitelevel_df = sitelevel_df[sitelevel_df['cov_score'] < thres].sort_values(by=['rname','pos'])
        bw = pyBigWig.open(f"{output_dir}/cov_score_not_pass_thres_{int(thres*10)}.bw", "wb")

        bw.addHeader(header,maxZooms=0)

        # index = sitelevel_df.index
        chroms = not_pass_sitelevel_df['rname'].values
        starts = np.array(not_pass_sitelevel_df['start'].values, dtype=np.int64)
        ends = np.array(not_pass_sitelevel_df['end'], dtype=np.int64)
        score = np.round(np.array(not_pass_sitelevel_df['cov_score'], dtype=np.float),4)
        if len(chroms.tolist()) > 0:
            bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=score.tolist(),validate=True)

        bw.close()
    # num reads

    bw = pyBigWig.open(f"{output_dir}/num_editted_reads.bw", "wb")
    bw.addHeader(header,maxZooms=0)
    chroms = sitelevel_df['rname'].values
    starts = np.array(sitelevel_df['start'].values, dtype=np.int64)
    ends = np.array(sitelevel_df['end'], dtype=np.int64)
    num_reads = np.array(sitelevel_df['methyl_cov'], dtype=np.float)
    bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=num_reads.tolist(),validate=True)
    bw.close()

    bw = pyBigWig.open(f"{output_dir}/num_all_reads.bw", "wb")
    bw.addHeader(header,maxZooms=0)
    chroms = sitelevel_df['rname'].values
    starts = np.array(sitelevel_df['start'].values, dtype=np.int64)
    ends = np.array(sitelevel_df['end'], dtype=np.int64)
    num_reads = np.array(sitelevel_df['cov'], dtype=np.float)
    bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=num_reads.tolist(),validate=True)
    bw.close()
        
        # bw = pyBigWig.open(f"{output_dir}/model_score_thres_{int(thres*10)}.bw", "wb")
        
        # bw.addHeader(header,maxZooms=0)
        # index = (sitelevel_df['model_score'] >= thres)
        # chroms = sitelevel_df.loc[index,'rname'].values
        # starts = sitelevel_df.loc[index,'start'].values
        # ends = np.array(sitelevel_df.loc[index,'end'], dtype=np.int64)
        # score = np.array(sitelevel_df.loc[index,'model_score'])

        # bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=score.tolist(),validate=True)

        # bw.close()
# def write_candidate_site(output_dir,candidate_path,reference_length_dict):
#     candidate_df = pd.read_csv(candidate_path,sep='\t',header=None)
#     candidate_df.columns = ['rname','pos','kmer','score','cov']
#     candidate_df['start'] = candidate_df['pos'] - 1
#     candidate_df['end'] = candidate_df['pos']
#     candidate_df = candidate_df.sort_values(['rname','pos'])
#     header = []
#     for rname in candidate_df['rname'].unique():
#         header.append((rname,reference_length_dict[rname]))
#     bw = pyBigWig.open(f"{output_dir}/candidate.bw", "wb")
#     bw.addHeader(header,maxZooms=0)
#     chroms = candidate_df.loc[:,'rname'].values
#     starts = candidate_df.loc[:,'start'].values
#     ends = np.array(candidate_df.loc[:,'end'], dtype=np.int64)
#     score = np.array(candidate_df.loc[:,'score'], dtype=np.float64)
#     bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=score.tolist(),validate=True)
#     bw.close()
def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--site_level_bed', type=str, default=None,required=False,help="Site level input bed")
    parser.add_argument('--mole_level_bed', type=str, default=None,required=False,help="Molecule level input bed")
    parser.add_argument('-ref','--reference_sequence', type=str, default=None,required=True,help="Reference sequence")
    parser.add_argument('--candidate_sites', type=str, default=None,required=False,help="Candidate site file path")
    parser.add_argument('--cDNA_mole_input', type=str, default=None,required=False,help="cDNA_mole_input")
    parser.add_argument('-o','--output_dir', type=str, default=None,required=True,help="Output path")
    parser.add_argument('-t','--threads', type=int, default=1,required=False,help="threads")
    args = parser.parse_args()
    reference_path = args.reference_sequence
    reference_fasta = pysam.FastaFile(reference_path)
    reference_id_dict = {}
    for i, name in zip(range(len(reference_fasta.references)), reference_fasta.references):
        reference_id_dict[name] = i
    reference_length_dict = {}
    for rname,length in zip(reference_fasta.references,reference_fasta.lengths):
        reference_length_dict[rname] = length
    output_dir = args.output_dir
    Path(output_dir).mkdir(exist_ok=True,parents=True)
    # use first 6 columns only
    if args.candidate_sites is not None:
        with open(args.candidate_sites,'r') as f:
            line = f.readline()
            if len(line) == 0:
                candidate_df = None
            else:
                candidate_df = pd.read_csv(args.candidate_sites,sep='\t',header=None)
                candidate_df.columns = ['rname','pos','kmer','score','cov']
                candidate_df['start'] = candidate_df['pos'] - 1
                candidate_df['end'] = candidate_df['pos']
                candidate_df = candidate_df.sort_values(['rname','pos'])
    else:
        candidate_df = None
#         write_candidate_site(output_dir,args.candidate_sites,reference_length_dict)
    if args.site_level_bed is not None:
        print('Start pre-compute site level visualization',flush=True)
        st_time = time.time()
        sitelevel_df = pd.read_csv(args.site_level_bed,sep='\t',header=None,usecols=range(6))
        sitelevel_df.columns = ['rname','pos','methyl_cov','cov','cov_score','model_score']
        sitelevel_df['start'] = sitelevel_df['pos'] - 1
        sitelevel_df['end'] = sitelevel_df['pos']
#         sitelevel_df[['rname','start','end','cov_score']].to_csv(f'{output_dir}/cov_score.bedgraph',sep='\t',header=None,index=None)
#         sitelevel_df[['rname','start','end','model_score']].to_csv(f'{output_dir}/model_score.bedgraph',sep='\t',header=None,index=None)
        # site level
        header = []
        for rname in sitelevel_df['rname'].unique():
            header.append((rname,reference_length_dict[rname]))
        pre_compute_site_level_thres(output_dir,header,sitelevel_df,candidate_df)
        print('DONE in '+str(time.time()-st_time)+' seconds',flush=True)
    if args.mole_level_bed is not None:
        print('Start pre-compute molecule level visualization',flush=True)
        st_time = time.time()
        # mole level
        mole_path = args.mole_level_bed
        mole_info_to_bam(mole_path,reference_id_dict,reference_path,output_dir,args.cDNA_mole_input,args.threads)
        print('DONE in '+str(time.time()-st_time)+' seconds',flush=True)
main()