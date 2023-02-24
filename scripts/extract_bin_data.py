import sys
import pysam
import numpy as np
import h5py
from collections import OrderedDict,deque
import re
import random
import argparse
from scipy.stats import skew, kurtosis
import os,glob,pickle
np.random.seed(1234)
import time
from pathlib import Path
import multiprocessing as mp

def debuginfoStr(info):
    print(info,flush=True)
    with open('/proc/self/status') as f:
        memusage = f.read().split('VmRSS:')[1].split('\n')[0][:-3]
    mem = int(memusage.strip())/1024
    pid = os.getpid()
    print(f'Current mem usage: {mem} MB reported by {pid}',flush=True)


def normaliza_signal(signal,summary_value):
    normalized = [(x-summary_value[0])/summary_value[1] for x in signal]
    return np.clip(normalized, -5, 5).tolist()

encoding_baseout = OrderedDict([
    ('A', 1),
    ('C', 2),
    ('G', 3),
    ('T',4),
    ('-', 0),
    ('N', 0),

])



compACGT = OrderedDict([('A','T'),('T','A')])

def reverse_compl(seq):
    res=''
    for s in seq:
        if s=='A':
           res+='T'
        if s=='T':
           res+='A'
        if s=='G':
           res+='C'
        if s=='C':
           res+='G'
        if s=='N':
           res+='N'
        if s=='-':
           res+='-'
    
    return res[len(res)::-1]

def frag_features(signal_frag_list,summary_value,featurenum=3):
    '''
    same feature as paper detection of DNA base modifications by deep recurrent neural network on Oxford Nanopore sequencing data
    '''
    featureMatrix = np.zeros([len(signal_frag_list),featurenum])
    if featurenum==3:
        for i in range(len(signal_frag_list)):
            if signal_frag_list[i] != '-':
               event = [float(x) for x in signal_frag_list[i].split(',')]
               event = normaliza_signal(event,summary_value)
               event_mean = np.mean(event)
               event_std = np.std(event)
               event_length = len(event)
            else:
               event_mean=0
               event_std=0
               event_length=0
            
            featureMatrix[i,0]=event_mean
            featureMatrix[i,1]=event_std
            featureMatrix[i,2]=event_length
    elif featurenum==5:
         for i in range(len(signal_frag_list)):
            if signal_frag_list[i] != '-':
               event = [float(x) for x in signal_frag_list[i].split(',')]
               event = normaliza_signal(event,summary_value)
               event_mean = np.mean(event)
               event_std = np.std(event)
               event_length = len(event)
               event_skew = skew(event)
               event_kurtosis = kurtosis(event)
            else:
               event_mean=0
               event_std=0
               event_length=0
               event_skew=0
               event_kurtosis=0
            
            featureMatrix[i,0]=event_mean
            featureMatrix[i,1]=event_std
            featureMatrix[i,2]=event_length
            featureMatrix[i,3]=event_skew
            featureMatrix[i,4]=event_kurtosis
    
    return featureMatrix

def extract_base_from_sam(molecule,ref_chr,sam_info):
    chrom=sam_info[molecule].reference_name
    sequence=sam_info[molecule].query_sequence
    st=sam_info[molecule].reference_start
    ed=sam_info[molecule].reference_end
    aligned=sam_info[molecule].aligned_pairs
    strand = sam_info[molecule].flag #0 (+) not 0 (-)
    seq_frags_call={}
    seq_frags_ref={}
    if sequence == None:
        return {},{}
    if coveragefile is not None:
        coveragepass_chr = {}
        page_number_st = st//page_size
        page_number_ed = ed//page_size
        list_of_files = []
        for page_number in range(page_number_st,page_number_ed+1):
            list_of_files += glob.glob(f'{temp_dir}/coveragepass_{chrom}_{page_number}_*.pkl')
        for file in list_of_files:
            with open(file,'rb') as f:
                coveragepass_chr.update(pickle.load(f))
    for i,j_ref in aligned: #i and j_ref are 0 based
        if j_ref is None:
            continue
        
        if candidate_file is not None:
            if candidateonly==1:
                #   print("candidate only\n")
                if chrom.split("|")[0]+"-"+str(j_ref+1) not in candidate_sites:
                    continue
            else:
                #  print("noncandidate only\n")
                if chrom.split("|")[0]+"-"+str(j_ref+1) in candidate_sites:
                    continue
        
        if coveragefile is not None:
            # print("use site coverage cutoff\n")
            if int(j_ref+1) not in coveragepass_chr:
                continue
        
        try:
            base = ref_chr[j_ref].upper() #aligned position start from 1
        except:
            #print("error!!!")
            #print(chrom)
            #print(j_ref+1)
            #print(len(ref_chr))
            print("Reference does not match bam file")
        
        if (base == center and strand ==0) or (base == compACGT[center] and strand !=0):
             center_index_ref = j_ref
             pos_left=np.max([center_index_ref-nt,0])
             pos_right=center_index_ref+nt
             refseq=ref_chr[pos_left:pos_right+1].upper()
             seq_frags_ref[center_index_ref+1] = refseq #can only be smaller than 11
             
             if len(refseq) > nt:
                  seq_frags_call[j_ref+1]=''
    
    
    ss={}
    for j in aligned:
        ss[j[1]]=j[0]
    
    #print(ss)
    #print(ed)
    for i in seq_frags_call: #i is 1-based not 0-based!
        region=[]
        for i1 in range(i-nt-1,i+nt):
            if i1<0 or i1>ed-1 or (i1 not in ss) or ss[i1]==None: #ss[il]=i in call
                region.append('-')
            elif i1<i+nt+1 and i1<ed-1 and type(ss[i1+1])==int:
                region.append(sequence[ss[i1]:ss[i1+1]])
            else:
                region.append(sequence[ss[i1]])
        #print(region)
        seq_frags_call[i]=''.join(region)	
    
    #print("debug seq_frags_call and ref")
    #print(seq_frags_call)
    #print(seq_frags_ref)
    return seq_frags_call,seq_frags_ref


def process_signal_for_one_read(read_id,mole,summary_one_read,read_information,ref_chr,sam_info):
    windowsize=2*nt+1
    maxoutputlen=windowsize+2 #if windowsize=11 this is 13
    center_index = int(windowsize/2) 
    # mole,seq_frags_call_one_read,seq_frags_ref,summary_one_read,read_information = args
    seq_frags_call_one_read,seq_frags_ref=extract_base_from_sam(read_id,ref_chr,sam_info)
    fragments,labels_ref,labels_call,info_list = [],[],[],[]
    if len(seq_frags_call_one_read)==0 or len(seq_frags_ref)==0:
          return fragments,labels_ref,labels_call,info_list
    
    labeltype,read_id,chrom = read_information['labeltype'],read_information['read_id'],read_information['chrom']
    
    for i in seq_frags_call_one_read:# here i is the position on ref
        if i not in seq_frags_ref:
             continue
        
        signal_frags=[]
        for j in range(i-nt-3,i+nt-2):
            if j not in mole:
                signal_frags.append('-')
            else:
                #a=Med(mole[str(j)],summary_one_read) #normalize raw signal and caculate std and mean for signle event not for window!
                signal_frags.append(','.join(mole[j])) #raw signal
        
        
        if signal_frags.count('-')>nt:
            #print("-\n")
            #count_delete+=1
            continue
        
        seq_frags_ref_ = seq_frags_ref[i]
        seq_frags_call_ = seq_frags_call_one_read[i]
        strand="+" #for + 0 for -
        if center_index>len(seq_frags_ref_):
              print("error:"+str(read_information['read_id'])+":"+str(read_information['chrom'])+"-"+str(i))
              #print(seq_frags_ref_)
        
        if seq_frags_ref_[center_index] == 'T':
                #print(seq_frags_ref_)
                #print(seq_frags_call_)
                seq_frags_ref_=reverse_compl(seq_frags_ref_)
                seq_frags_call_=reverse_compl(seq_frags_call_)
                signal_frags.reverse()
                strand='-'
        
        if candidate_file is not None and candidateonly==1:
             info_ = "\t".join([labeltype+":"+str(candidate_sites[chrom.split("|")[0]+"-"+str(i)]),read_id,chrom,str(i),strand])
        else:
             info_ = "\t".join([labeltype,read_id,chrom,str(i),strand])
        
        
        features = frag_features(signal_frags,summary_one_read,FEATURENUM)
        #if "N" in seq_frags_ref_:
        #     print(seq_frags_ref_)
        #if len(seq_frags_ref_)>maxoutputlen:
        #   print("error at readid="+read_id+str(chrom)+" "+str(center_index))
        #
        window_base_ref = [encoding_baseout[b] for b in seq_frags_ref_]
        window_base_call = [encoding_baseout[b] for b in seq_frags_call_]
        if len(window_base_call)>maxoutputlen:
            #kept everything
            window_base_call=window_base_call[:maxoutputlen]
        else:
            window_base_call=window_base_call + [0] * (maxoutputlen - len(window_base_call))
        
        if len(window_base_ref)> maxoutputlen:
             window_base_ref=window_base_ref[:maxoutputlen]
        else:
              window_base_ref=window_base_ref + [0] * (maxoutputlen - len(window_base_ref)) #0 is for -
        
        fragments.append(features)
        labels_ref.append(window_base_ref)
        labels_call.append(window_base_call)
        info_list.append(info_)
    
    return fragments,labels_ref,labels_call,info_list

def get_summary_and_sam_info(read_id_set):
    sam_info={} #load all sam_info into memory
    with pysam.AlignmentFile(samIn, "rb") as samfile:
        for line in samfile:
            molecule=line.query_name
            if molecule not in read_id_set:
                continue
            if molecule not in sam_info:
                 sam_info[molecule]=line
            else:
                 print("not unique mappint for "+molecule+" !!!\n")
                 exit(0) 
                 #sam_info[molecule].append(line)
    
    #print("load mapping information over\n")
    #debuginfoStr("load mapping information over")
    #G1+G2=369+352
    #get guppy summary info
    summary={}
    with open(summaryIn) as fidIn1:
        names=next(fidIn1).strip().split('\t')
        read_id=names.index('read_id')
        med_id=names.index('median_template')
        mad_id=names.index('mad_template')
        
        for line in fidIn1:
            if "filename" in line:
                continue
            info=line.strip().split('\t')
            
            if info[read_id] not in read_id_set:
                continue
            
            summary[info[read_id]]=[float(info[med_id]),float(info[mad_id])]
    
    #debuginfoStr("load summary")
    return sam_info,summary

def process_signals(args):
    os.nice(10)
    msgs_queue,eventIn,start_file_pos,start_read_index,num_reads_assigned,ref = args
    mole_buffer,read_information_buffer = [],[]
    read_id_set = set(all_reads_id[start_read_index:start_read_index+num_reads_assigned])
    all_reads_id.clear()
    sam_info,summary = get_summary_and_sam_info(read_id_set)
    
    with open(eventIn) as fidIn5:
        fragments,labels_ref,labels_call,info_list = [],[],[],[]
        header = next(fidIn5).strip().split("\t")
        position_index = header.index("position")
        contig_index = header.index("contig")
        read_name_index = header.index("read_name")
        raw_signal_index = header.index("samples")
        event_index = header.index("event_index")
        mole={} #pos list of signals
        read_id=''#readid
        chrom=''
        pre_event_index=0
        num_reads_processed = 0
        fidIn5.seek(start_file_pos)
        for line in fidIn5:
            if (num_reads_processed >= num_reads_assigned):
                break
            info=line.strip().split('\t')
            if len(info)<16: #make sure pipline can continue.
                continue
            
            raw=info[raw_signal_index] #raw signal -1
            next_read_id=info[read_name_index] #readid 3 
            pos = int(info[position_index]) #info[1] is position on ref
            next_event_index =int(info[event_index])
            if read_id!=next_read_id:
                if mole!={}:
                    mole_buffer.append(mole)
                    read_information_buffer.append({'labeltype':labeltype,'read_id':read_id,'chrom':chrom})
                    num_reads_processed += 1
                    if ((len(mole_buffer) + 1) > MAX_LEN_MOLE_BUFFER) or (num_reads_processed >= num_reads_assigned):
                        fragments,labels_ref,labels_call,info_list = [],[],[],[]
                        read_ids = [read_information['read_id'] for read_information in read_information_buffer]
                        assert len(mole_buffer) == len(read_information_buffer)
                        for read_id,mole,read_information in zip(read_ids,mole_buffer,read_information_buffer):
                            chrom = read_information['chrom']
                            if chrom not in ref:
                                with open(f'{temp_dir}/ref_{chrom}.txt','r') as f:
                                    ref[chrom] = f.read()
                            
                            
                            
                            
                            
                            res = process_signal_for_one_read(read_id,mole,summary[read_id],read_information,ref[chrom],sam_info)
                            if res[0] !=[]:
                                fragments += res[0]
                                labels_ref += res[1]
                                labels_call += res[2]
                                info_list += res[3]
                        msgs_queue.put((fragments,labels_ref,labels_call,info_list,len(read_ids)))
                        for lst in [mole_buffer,read_information_buffer,fragments,labels_ref,labels_call,info_list]:
                            lst.clear()
                        
                read_id=next_read_id
                mole={}
                chrom=info[contig_index]
            else:
                if pos not in mole:
                    mole[pos]=deque([])
                if next_event_index > pre_event_index: #strand ==0 
                        mole[pos].append(raw)
                else: #strand !=0
                        mole[pos].appendleft(raw)
                
                pre_event_index=next_event_index
def listener(msgs_queue,cachefile,windowsize,featuredim,maxoutputlen):
    h5f = h5py.File(cachefile, 'w')
    h5f.create_dataset('X', (0,windowsize,featuredim),maxshape=(None,windowsize,featuredim), dtype='float32')#features has 3+4 
    h5f.create_dataset('y_ref', (0,maxoutputlen),maxshape=(None,maxoutputlen), dtype='float32')
    h5f.create_dataset('y_call', (0,maxoutputlen),maxshape=(None,maxoutputlen), dtype='float32')
    h5f.create_dataset('info', (0,),maxshape=(None,), dtype=h5py.special_dtype(vlen=str))#features has 3+4 
    saved_size = 0
    # for debuging
    total_num_reads,num_reported = 0,0
    
    while True:
        msg = msgs_queue.get()
        if msg == 'kill':
            break
        else:
            fragments,labels_ref,labels_call,info_list,num_reads = msg
            # for debuging
            total_num_reads += num_reads
            #if total_num_reads // 1e4 > num_reported:
            #    print(f'Processed {total_num_reads} reads',flush=True)
            #    print(f'Processed {saved_size} signals',flush=True)
            #    num_reported += 1
            
            h5f['info'].resize((saved_size+len(info_list),))
            h5f['info'][saved_size:,]=info_list
            h5f['X'].resize((saved_size+len(fragments),windowsize,featuredim))
            h5f['X'][saved_size:,:,:]=np.asarray(fragments)
            h5f['y_ref'].resize((saved_size+len(labels_ref),maxoutputlen))
            h5f['y_ref'][saved_size:,:]=labels_ref
            h5f['y_call'].resize((saved_size+len(labels_call),maxoutputlen))
            h5f['y_call'][saved_size:,:]=labels_call
            h5f.flush()
            saved_size+=len(fragments) 
##def savetohdf5()
    #python3 extract_bin_data.py "A" 5 "reference.fa" "test_summary.txt" "test_samefile.sam" "test_eventfile.txt" "I" "test_output.hdf5"
def dump_ref(Ref,dump_to_disk):
    ref_dict = {}
    with open(Ref) as fid1:
        seq=''
        chrom=next(fid1).strip().split()[0][1:]
        for line in fid1:
            #print(line)
            if line.startswith('>'):
                if dump_to_disk:
                    with open(f'{temp_dir}/ref_{chrom}.txt','w') as f:
                        f.write(seq)
                else:
                    ref_dict[chrom] = seq
                chrom=line.strip().split()[0][1:]
                seq=''
            else:
                seq+=line.strip()
        if dump_to_disk:
            with open(f'{temp_dir}/ref_{chrom}.txt','w') as f:
                f.write(seq)
        else:
            ref_dict[chrom] = seq
    return ref_dict


def dump_coveragepass_iteration(args):
    coveragefile,coveragecutoff,start_file_pos,num_lines_per_process,process_id = args
    coveragepass_chr = {}
    curr_chr = ''
    curr_page_number = -1
    with open(coveragefile) as file:
        file.seek(start_file_pos)
        line_num_ct = 0
        for line in file:
            if line_num_ct >= num_lines_per_process:
                break
            chr_ = line.split()[0]
            pos_ = int(line.split()[1])
            coverage=float(line.split()[-1])
            line_num_ct += 1
            if coverage>=0:
                page_number = pos_ // page_size
                if ((page_number != curr_page_number) and (curr_page_number != -1)) or ((chr_ != curr_chr) and (curr_chr != '')):
                    with open(f'{temp_dir}/coveragepass_{curr_chr}_{curr_page_number}_{process_id}.pkl','wb') as f:
                        pickle.dump(coveragepass_chr,f)
                    coveragepass_chr = {}
                coveragepass_chr[pos_] = True
                curr_page_number = page_number
                curr_chr = chr_
    with open(f'{temp_dir}/coveragepass_{curr_chr}_{curr_page_number}_{process_id}.pkl','wb') as f:
        pickle.dump(coveragepass_chr,f)
def get_coveragepass_args(coveragefile,coveragecutoff,threads):
    line_offset = []
    offset = 0
    with open(coveragefile) as file:
        for line in file:
            line_offset.append(offset)
            offset += len(line)
    
    num_lines_per_process, extra = divmod(len(line_offset), threads)
    if extra:
        num_lines_per_process += 1
    list_of_args = []
    for i in range(threads):
        list_of_args.append((coveragefile,coveragecutoff,line_offset[i*num_lines_per_process],num_lines_per_process,i))
    return list_of_args

def dump_coveragepass(coveragefile,coveragecutoff,threads):
    list_of_args = get_coveragepass_args(coveragefile,coveragecutoff,threads)
    pool = mp.Pool(threads)
    futures = []
    for args in list_of_args:
        futures.append(pool.apply_async(dump_coveragepass_iteration,(args,)))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
    #print('coveragepass dumped',flush=True)
def get_eventIn_offset(eventIn,threads):
    all_reads_id = []
    with open(eventIn) as fidIn5:
        line_offset = []
        line = fidIn5.readline()
        offset = fidIn5.tell()
        header = line.strip().split("\t")
        read_name_index = header.index("read_name")
        read_id=''
        for line in fidIn5:
            info=line.strip().split('\t')
            offset += len(line)
            if len(info)<16: #make sure pipline can continue.
                continue
            next_read_id=info[read_name_index]
            if read_id != next_read_id:
                line_offset.append(offset)
                read_id = next_read_id
                all_reads_id.append(read_id)
    num_reads = len(line_offset)
    assert len(all_reads_id) == num_reads
    #print(f'Total number of reads {num_reads}',flush=True)
    chunksize, extra = divmod(num_reads, threads)
    if extra:
        chunksize += 1
    return [(line_offset[i*chunksize],i*chunksize) for i in range(threads)],all_reads_id,chunksize
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--center', type=str, default="A",
                        help="center")
    parser.add_argument('--nt', type=int, default=5,
                        help="number of neighbors")
    parser.add_argument('--Ref', type=str, default=None,
                        help="reference genome")
    parser.add_argument('--summaryIn', type=str, default=None,
                        help="summary file")
    parser.add_argument('--samIn', type=str, default=None,
                        help="samfile file")
    parser.add_argument('--eventIn', type=str, default=None,
                        help="event file")
    parser.add_argument('--labeltype', type=str, default=None,
                        help="labeltype")
    parser.add_argument('--cachefile', type=str, default=None,
                        help="cachefile")
    parser.add_argument('--THREADS', type=int, default=6,
                        help="number of threads")
    
    parser.add_argument('--maxbuffer_size', type=int, default=1000,
                        help="maxbuffer_size")
    
    parser.add_argument('--candidate_file', type=str, default=None,
                        help="candidate site file")
    
    parser.add_argument('--featurenum', type=int, default=3,
                        help="featurenum. 3: mean, std, length")
    
    parser.add_argument('--candidateonly', type=int, default=1,
                        help="candidateonly. 1: candidateonly otherwise noncandidate (candite_file must not None)")
    
    parser.add_argument('--coveragefile', type=str, default=None,
                        help="coveragefile, only site in this file can be passed")
    
    parser.add_argument('--coveragecutoff', type=int, default=10,
                        help="site coverage threshold")
    parser.add_argument('--ref_dump_position', type=str, default='memory',
                        help="Dump reference files in disk or memory")
    parser.add_argument('--temp_dir', type=str, default=None,
                        help="temp folder")
    
    args = parser.parse_args()
    #for k, v in vars(args).items():
    #    
    #    #print(k, ':', v)
    
    center=args.center
    nt=args.nt
    Ref=args.Ref
    summaryIn=args.summaryIn
    samIn=args.samIn
    eventIn=args.eventIn
    labeltype=args.labeltype
    cachefile=args.cachefile
    THREADS=args.THREADS
    candidate_file = args.candidate_file 
    MAX_LEN_MOLE_BUFFER = args.maxbuffer_size
    FEATURENUM = args.featurenum
    candidateonly=args.candidateonly
    if args.ref_dump_position == 'disk':
        dump_to_disk = True
    else:
        dump_to_disk = False
    temp_dir = args.temp_dir
    if temp_dir is None:
         temp_dir = os.path.join(os.path.dirname(args.cachefile),os.path.basename(args.cachefile)+"temp")
    
    start_time = time.time()
    Path(temp_dir).mkdir(parents=True,exist_ok=True)
    for file in glob.glob(f'{temp_dir}/*'):
        Path(file).unlink()
    
    featuredim=FEATURENUM
    windowsize=2*nt+1
    maxoutputlen=windowsize+2 #if windowsize=11 this is 13
    center_index = int(windowsize/2) 
    coveragefile=args.coveragefile
    page_size = 16384
    if coveragefile is not None:
        dump_coveragepass(coveragefile,args.coveragecutoff,THREADS)
        #debuginfoStr("dump coveragepass file")
    ref_dict = dump_ref(Ref,dump_to_disk)
    #debuginfoStr("read reference over")
    if candidate_file is not None:
        candidate_sites = {}
        with open(candidate_file) as file:
           for line in file:
               chr_ = line.split()[0]
               pos_ = line.split()[1]
               if len(line.split("\t"))<9:
                  candidate_sites[chr_+"-"+pos_]=float(line.split("\t")[3])
               else:
                  candidate_sites[chr_+"-"+pos_]=float(line.split("\t")[8])
    
    #debuginfoStr("load candidate_file")
    list_of_offsets,all_reads_id,chunksize = get_eventIn_offset(eventIn,THREADS)
    #debuginfoStr("get eventIn offset")
    manager = mp.Manager()
    pool = mp.Pool(THREADS+1)
    msgs_queue = manager.Queue() 
    list_of_args = []
    for offset,start_read_index in list_of_offsets:
        list_of_args.append((msgs_queue,eventIn,offset,start_read_index,chunksize,ref_dict))
    
    watcher = pool.apply_async(listener, args=(msgs_queue,cachefile,windowsize,featuredim,maxoutputlen))
    futures = []
    for args in list_of_args:
        futures.append(pool.apply_async(process_signals,(args,)))
    for future in futures:
        future.get()
    msgs_queue.put('kill')
    watcher.get()
    pool.close()
    pool.join()
    #count_delete = 0
    for file in glob.glob(f'{temp_dir}/*'):
        Path(file).unlink()
    Path(f'{temp_dir}').rmdir()
print("--- succucessful---")
print("--- %s seconds ---" % (time.time() - start_time))
#python3 extract_bin_data.py "A" 5 "reference.fa" "test_summary.txt" "test_samefile.sam" "test_eventfile.txt" "I" "test_output.hdf5"