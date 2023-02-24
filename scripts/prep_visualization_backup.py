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
import array
import re
import bisect
def get_color(score):
    score = float(score)
    r_val = int((1-score) * 255)
    b_val = int(score * 255)
    return f"{r_val},0,{b_val}"
# def get_cigar(sites_list,score_list,thres):
#     cigar_list = []
#     reference_start = sites_list[0]
#     reference_start_score = score_list[0]
#     previous_pos = reference_start
#     if reference_start_score >= thres:
#         cigar_list.append((0,1))
#     else:
#         cigar_list.append((2,1))
#     if len(sites_list) > 1:
#         for pos,score in zip(sites_list[1:],score_list[1:]):
#             distance = pos - previous_pos + 1 - 2
#             if distance > 0:
#                 cigar_list.append((3,distance))
#                 if score >= thres:
#                     cigar_list.append((0,1))
#                 else:
#                     cigar_list.append((2,1))
#             elif distance == 0:
#                 if score >= thres:
#                     cigar_list.append((0,1))
#                 else:
#                     cigar_list.append((2,1))
#             previous_pos = pos
#     return cigar_list
def get_cigar(sites_list,score_list,thres):
    cigar_list = []
    reference_start = sites_list[0]
    reference_start_score = score_list[0]
    previous_pos = reference_start
    if reference_start_score >= thres:
        cigar_list.append((7,1))
    else:
        cigar_list.append((7,1))
    if len(sites_list) > 1:
        for pos,score in zip(sites_list[1:],score_list[1:]):
            distance = pos - previous_pos + 1 - 2
            if distance > 0:
                cigar_list.append((3,distance))
                if score >= thres:
                    cigar_list.append((7,1))
                else:
                    cigar_list.append((7,1))
            elif distance == 0:
                if score >= thres:
                    cigar_list.append((7,1))
                else:
                    cigar_list.append((7,1))
            previous_pos = pos
    return cigar_list
def get_query_qualities(score_list):
    query_qualities = []
    color_thres = 0.5
    minQ = 5
    maxQ = 20
    color_thres = [i/10 for i in range(1,11)]
    rgba_list = [[0, 0, 255, 1.0],
     [0, 0, 255, 0.5],
     [0, 0, 255, 0.3],
     [255, 0, 0, 0.4],
     [255, 0, 0, 0.5],
     [255, 0, 0, 0.7],
     [255, 0, 0, 0.9],
     [195, 1, 1, 1.0],
     [148, 0, 0, 0.9],
     [148, 0, 0, 1.0]]
    for score in score_list:
        if score < 0:
            score = 0
        if score > 1:
            score = 1
        alpha = rgba_list[bisect.bisect_left(color_thres,score)][3]
        Q = int(minQ+(maxQ-minQ)*alpha)
        query_qualities.append(Q)
    return array.array('B',query_qualities)
def get_query_sequence(score_list):
    query_sequences = []
    color_thres = [i/10 for i in range(1,11)]
    base_list = ['A' for i in range(3)] + ['T' for i in range(4)] +['C'] + ['G' for i in range(2)]
    for score in score_list:
        if score < 0:
            score = 0
        if score > 1:
            score = 1
        base = base_list[bisect.bisect_left(color_thres,score)]
        query_sequences.append(base)

    return ''.join(query_sequences)
def get_line_marker(mole_path,threads):
    '''
    Split the file into THREADS chunks
    !Split by read
    '''
    file_stats = os.stat(mole_path)
    total_bytes = file_stats.st_size
    chunksize, extra = divmod(total_bytes, threads)
    if extra:
        chunksize += 1
    byte_marker = []
    previous_read_name = ''
    with open(mole_path,'r') as f:
        for i in range(threads):
            if i == 0:
                start_offset = 0
                for line in f:
                    if line[0] != '#':
                        break
                    start_offset += len(line)
            else:
                f.seek(i*chunksize)
                f.readline()
                while True:
                    line = f.readline()
                    read_name = line.split('\t')[1]
                    if previous_read_name == '':
                        # first (possible) incomplete line. Just read it.
                        previous_read_name = read_name
                        start_offset = f.tell()
                    else:
                        if previous_read_name != read_name:
                            # new read
                            previous_read_name = ''
                            break
                        else:
                            # still current read
                            previous_read_name = read_name
                            start_offset = f.tell()
            byte_marker.append(start_offset)
    byte_marker.append(total_bytes)
    line_marker = []
    for i in range(threads):
         line_marker.append((byte_marker[i],byte_marker[i+1]))
    return line_marker
def write_to_bam(old_read_id,old_read_sites,old_read_info,reference_id_dict,thres,reference_type,output_handle):
    score_list = []
    sites_list = []
    old_read_sites = sorted(old_read_sites,key=lambda x:x[0])
    for pos,score in old_read_sites:
        score_list.append(float(score))
        sites_list.append(pos)
    seg = pysam.AlignedSegment()
    seg.query_name = str(old_read_id)
    seg.flag = 0 if old_read_info['strand'] == '+' else 16
    seg.cigar = get_cigar(sites_list,score_list,thres)
    seg.reference_id = reference_id_dict[old_read_info['chr_name']]
    seg.reference_start = sites_list[0] - 1
    seg.mapping_quality = 60
    seg.query_sequence = get_query_sequence(score_list)
    seg.query_qualities = get_query_qualities(score_list)
    score_str = str(round(np.mean(score_list),2))
#                 seg.tags = [('SC',get_color(np.mean(score_list)),'Z'),('ML',array.array('f',score_list),'f')]
    tags = [('SC',get_color(np.mean(score_list)),'Z'),('ML',f'Score:{score_str}','Z')]
    if reference_type == 'transcriptome':
        if old_read_info['isoform'] != '':
            tags.append(('TT',old_read_info['isoform'],'Z'))
    seg.tags = tags
    output_handle.write(seg)
def annotation_to_bam(gtf_path,reference_path,temp_folder):
    transcript_dict = {}
    with open(gtf_path) as f:
        for line in f:
            fields = line.split('\t')
            if fields[2] == 'exon':
                gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
                isoform_name = re.findall('transcript_id "([^"]*)"', fields[8])[0]
                start_pos = int(fields[3])
                end_pos = int(fields[4])
                if isoform_name not in transcript_dict:
                    transcript_dict[isoform_name]  = {'start_pos':[],'end_pos':[]}
                transcript_dict[isoform_name]['start_pos'].append(start_pos)
                transcript_dict[isoform_name]['end_pos'].append(end_pos)
                transcript_dict[isoform_name]['chr'] = fields[0]
                transcript_dict[isoform_name]['strand'] = fields[6]
    reference_fasta = pysam.FastaFile(reference_path)
    annot_bam_path = f'{temp_folder}/annotation.bam'
    fout = pysam.AlignmentFile(annot_bam_path,'wb',reference_names=reference_fasta.references,reference_lengths=reference_fasta.lengths,threads=1)
    for isoform in transcript_dict:
        start = min(transcript_dict[isoform]['start_pos'])-1
        end = max(transcript_dict[isoform]['end_pos'])
        length = end - start + 1
        seg = pysam.AlignedSegment()
        seg.query_name = isoform
        seg.flag = 0 if transcript_dict[isoform]['strand']=='+' else 16
        seg.cigar = [(2,length)]
        seg.reference_id = reference_fasta.references.index(transcript_dict[isoform]['chr'])
        seg.reference_start = start
        seg.mapping_quality = 60
        tags = []
        tags.append(('TT',isoform,'Z'))
        seg.tags = tags
        fout.write(seg)
    fout.close()
    return annot_bam_path
def mole_info_to_bam_worker(mole_path,reference_type,reference_id_dict,reference_path,mole_output_dir,worker_id,thres,marker):
    start_file_pos,end_file_pos = marker
    # read_isoform_mapping = {}
    # if reference_type is not None:
    #     with open(reference_type) as f:
    #         for line in f:
    #             fields = line.split('\t')
    #             read_id = fields[1]
    #             isoform = fields[2]
    #             read_isoform_mapping[read_id] = isoform
    reference_fasta = pysam.FastaFile(reference_path)
    num_processsed_reads = 0
    old_read_id = ''
    old_read_sites = []
    old_read_info = {}
    temp_folder = f'{mole_output_dir}/temp/'
    Path(temp_folder).mkdir(exist_ok=True,parents=True)
    output_handle = pysam.AlignmentFile(f'{temp_folder}/{worker_id}.bam', "wb",reference_names=reference_fasta.references,reference_lengths=reference_fasta.lengths,threads=1)
    with open(mole_path,'r') as f:
        f.seek(start_file_pos)
        while True:
            line = f.readline()
            if not line:
                break
            if f.tell() > end_file_pos:
                break
            fields = line.split('\t')
            if reference_type == 'transcriptome':
                read_id,isoform,strand,score = fields[1],fields[2],fields[4],fields[5]
                if len(fields) == 8:
                    chr_name,start_pos = fields[6],fields[7]
                else:
                    continue
            elif reference_type == 'genome':
                [_,read_id,chr_name,start_pos,strand,score] = fields
#             if float(score) < thres:
#                 continue
            if chr_name not in reference_id_dict:
                continue
            start_pos = int(start_pos)
            if read_id != old_read_id and old_read_id != '':
                num_processsed_reads += 1
                write_to_bam(old_read_id,old_read_sites,old_read_info,reference_id_dict,thres,reference_type,output_handle)
                # new read
                old_read_id = read_id
                old_read_sites = []
                old_read_info = {}
            old_read_id = read_id
            old_read_sites.append((start_pos,score))   
            old_read_info = {'read_id':read_id,'chr_name':chr_name,'strand':strand}
            if reference_type == 'transcriptome':
                old_read_info['isoform'] = isoform
        # if reach last line of file
        if len(old_read_info) != 0:
            write_to_bam(old_read_id,old_read_sites,old_read_info,reference_id_dict,thres,reference_type,output_handle)
            output_handle.close()
def mole_info_to_bam(mole_path,reference_id_dict,reference_path,gtf_path,mole_output_dir,reference_type,threads):
    line_marker = get_line_marker(mole_path,threads)
    temp_folder = f'{mole_output_dir}/temp/'
    thres = 0.5
    pool = mp.Pool(threads+1)
    futures = []
    all_bam_files = []
    for marker,i in zip(line_marker,range(threads)):
        all_bam_files.append(f'{temp_folder}/{i+1}.bam')
        futures.append(pool.apply_async(mole_info_to_bam_worker,(mole_path,reference_type,reference_id_dict,reference_path,mole_output_dir,i+1,thres,marker,)))
    annot_bam_path = annotation_to_bam(gtf_path,reference_path,temp_folder)
    all_bam_files.append(annot_bam_path)
    for future in futures:
        future.get()
    pool.close()
    pool.join()
    sorted_bam = f'{mole_output_dir}/molecule_level.sorted.bam'
    bam = f'{mole_output_dir}/molecule_level.bam'
    params = ['-@',str(threads),'-f',bam]+all_bam_files
    pysam.merge(*params)
    pysam.sort('-@',str(threads),"-o",sorted_bam,bam)
    # pysam.sort('-@',str(threads),'-t',"TT","-n","-o",sorted_bam,bam)
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
            if len(chroms.tolist()) != 0:
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
            if len(chroms.tolist()) != 0:
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
            if len(chroms.tolist()) != 0:
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
    if len(chroms.tolist()) != 0:
        bw.addEntries(chroms.tolist(), starts.tolist(), ends=ends.tolist(), values=num_reads.tolist(),validate=True)
    bw.close()

    bw = pyBigWig.open(f"{output_dir}/num_all_reads.bw", "wb")
    bw.addHeader(header,maxZooms=0)
    chroms = sitelevel_df['rname'].values
    starts = np.array(sitelevel_df['start'].values, dtype=np.int64)
    ends = np.array(sitelevel_df['end'], dtype=np.int64)
    num_reads = np.array(sitelevel_df['cov'], dtype=np.float)
    if len(chroms.tolist()) != 0:
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
def isoform_site_to_bam(gtf_path,reference_path,trans_site_path,output_dir,threads):
    bam = f'{output_dir}/isoform_site.bam'
    sorted_bam = f'{output_dir}/isoform_site.sorted.bam'
    reference_fasta = pysam.FastaFile(reference_path)
    gene_dict = {}
    with open(gtf_path) as f:
        for line in f:
            fields = line.split('\t')
            if fields[2] == 'exon':
                gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
                start_pos = int(fields[3])
                end_pos = int(fields[4])
                if gene_name not in gene_dict:
                    gene_dict[gene_name]  = {'start_pos':[],'end_pos':[]}
                gene_dict[gene_name]['start_pos'].append(start_pos)
                gene_dict[gene_name]['end_pos'].append(end_pos)
                gene_dict[gene_name]['chr'] = fields[0]
                gene_dict[gene_name]['strand'] = fields[6]
    fout = pysam.AlignmentFile(bam,'wb',reference_names=reference_fasta.references,reference_lengths=reference_fasta.lengths,threads=1)
    for gene in gene_dict:
        start = min(gene_dict[gene]['start_pos'])-1
        end = max(gene_dict[gene]['end_pos'])
        length = end - start + 1
        seg = pysam.AlignedSegment()
        seg.query_name = gene
        seg.flag = 0
        seg.cigar = [(2,length)]
        seg.reference_id = reference_fasta.references.index(gene_dict[gene]['chr'])
        seg.reference_start = start
        seg.mapping_quality = 60
        tags = []
        tags.append(('GE',gene,'Z'))
        seg.tags = tags
        fout.write(seg)
    df = pd.read_csv(trans_site_path,sep='\t',header=None,usecols=[0,4,5,7,8,9])
    df.columns = ['isoform','cov_score','model_score','rname','pos','gene']
    df.groupby('isoform').apply(lambda rows:rows_to_read(rows,fout,reference_fasta))
    fout.close()
    pysam.sort('-@',str(threads),"-o",sorted_bam,bam)
    pysam.index(sorted_bam)
def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--site_level_bed', type=str, default=None,required=False,help="Site level input bed")
    parser.add_argument('--trans_site_level_bed', type=str, default=None,required=False,help="Isoform site level input bed")
    parser.add_argument('--mole_level_bed', type=str, default=None,required=False,help="Molecule level input bed")
    parser.add_argument('-ref','--reference_sequence', type=str, default=None,required=True,help="Reference sequence")
    parser.add_argument('--candidate_sites', type=str, default=None,required=False,help="Candidate site file path")
    parser.add_argument('--annotation_gtf', type=str, default=None,required=False,help="Annotation_gtf_path")
    parser.add_argument('--reference_type', type=str, default='transcriptome',required=False,help="The type of reference used [transcriptome,genome] [default:transcriptome]")
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
    if args.mole_level_bed is not None:
        print('Start pre-compute molecule level visualization',flush=True)
        st_time = time.time()
        # mole level
        mole_path = args.mole_level_bed
        mole_info_to_bam(mole_path,reference_id_dict,reference_path,args.annotation_gtf,output_dir,args.reference_type,args.threads)
        print('DONE in '+str(time.time()-st_time)+' seconds',flush=True)
    if args.site_level_bed is not None:
        print('Start pre-compute site level visualization',flush=True)
        st_time = time.time()
        sitelevel_df = pd.read_csv(args.site_level_bed,sep='\t',header=None,usecols=range(6))
        sitelevel_df.columns = ['rname','pos','methyl_cov','cov','cov_score','model_score']
        sitelevel_df = sitelevel_df[sitelevel_df['rname'].isin(set(reference_length_dict.keys()))]
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
    if args.trans_site_level_bed is not None:
        print('Start pre-compute isoform specific site level visualization',flush=True)
        isoform_site_to_bam(args.annotation_gtf,reference_path,args.trans_site_level_bed,output_dir,args.threads)
main()