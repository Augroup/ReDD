## GPU



#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1 --cpus-per-task=6 --gpus=1 --mem=60GB
#SBATCH --account=kinfai0 --partition=spgpu

source activate /scratch/kinfai_root/kinfai0/haorli/software/ReDD_env
cd /scratch/kinfai_root/kinfai0/haorli/ReDD_rna004/HEK293T
/gpfs/accounts/kinfai_root/kinfai0/haorli/software/dorado-1.0.2-linux-x64/bin/dorado basecaller hac /scratch/kinfai_root/kinfai0/haorli/KaiWang_dRNA/HEK293T_directRNA/pod5_pass/ --emit-fastq > pod5_pass.fastq
blue-crab p2s /scratch/kinfai_root/kinfai0/haorli/KaiWang_dRNA/HEK293T_directRNA/pod5_pass/ -o pod5_pass.blow5
/scratch/kinfai_root/kinfai0/haorli/software/f5c_ReDD/f5c index pod5_pass.fastq --slow5 pod5_pass.blow5
minimap2 -ax splice -uf -k14 --secondary=no -t6 genome.fasta pod5_pass.fastq | samtools sort -@ 6 - -o pod5_pass.bam && samtools index pod5_pass.bam
/scratch/kinfai_root/kinfai0/haorli/software/f5c_ReDD/f5c eventalign --rna --min-mapq 0 --pore rna004 -r pod5_pass.fastq -b pod5_pass.bam -g /nfs/turbo/umms-kinfai/yunhwang/Reference/Human/Genome/hg38/genome.fa --slow5 pod5_pass.blow5 \
-o window_size_9_hg38 \
--redd-candidate /nfs/turbo/umms-kinfai/haorli/20240314_ReDD_result_data/figure2a/reditools2_candidates/HEK293T_WT.candidate_sites.tab \
--redd --redd-window-size 9 --scale-events -B 20M -K 1024





/scratch/kinfai_root/kinfai0/haorli/software/f5c_ReDD/f5c eventalign --rna --min-mapq 0 --pore rna004 -r pod5_pass.fastq -b pod5_pass.bam -g genome.fasta --slow5 pod5_pass.blow5 \
-o pod5_pass_f5c_CPU \
--redd-candidate /nfs/turbo/umms-kinfai/haorli/20240314_ReDD_result_data/figure2a/reditools2_candidates/HEK293T_WT.candidate_sites.tab \
--print-read-names --redd --redd-window-size 11 --scale-events --disable-cuda=yes


run eventalign --rna --min-mapq 0 --pore rna004 -r pod5_pass.fastq.bgz -b pod5_pass.bam -g /scratch/kinfai_root/kinfai0/haorli/ReDD_rna004/reference/genome.fa.bgz --slow5 pod5_pass.blow5 -o window_size_9_hg38 --redd-candidate /nfs/turbo/umms-kinfai/haorli/20240314_ReDD_result_data/figure2a/reditools2_candidates/HEK293T_WT.candidate_sites.tab --redd --redd-window-size 9 --scale-events -B 20M -K 10 --debug-break 1
python /scratch/kinfai_root/kinfai0/haorli/software/ReDD/scripts/predict.py --device 'cuda' --threads 12 --batch_size 10000 --weights_dir /scratch/kinfai_root/kinfai0/haorli/software/ReDD/scripts/models/HEK293T_specific.ckpt --test_cachedata window_size_9_hg38.noncandidate.hdf5 --test_outputfile ./window_size_9_hg38.noncandidate.prediction.raw.txt
python /scratch/kinfai_root/kinfai0/haorli/software/ReDD/scripts/predict.py --device 'cuda' --threads 12 --batch_size 10000 --weights_dir /scratch/kinfai_root/kinfai0/haorli/software/ReDD/scripts/models/HEK293T_specific.ckpt --test_cachedata window_size_9_hg38.candidate.hdf5 --test_outputfile ./window_size_9_hg38.candidate.prediction.raw.txt



#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1 --cpus-per-task=6 --gpus=1 --mem=150GB
#SBATCH --account=kinfai0 --partition=spgpu

source activate /scratch/kinfai_root/kinfai0/haorli/software/ReDD_env
cd /scratch/kinfai_root/kinfai0/haorli/ReDD_rna004/K562
/gpfs/accounts/kinfai_root/kinfai0/haorli/software/dorado-1.0.2-linux-x64/bin/dorado basecaller hac /scratch/kinfai_root/kinfai0/haorli/KaiWang_dRNA/K562_directRNA/pod5_pass/ --emit-fastq > pod5_pass.fastq
bgzip < pod5_pass.fastq > pod5_pass.fastq.bgz
blue-crab p2s /scratch/kinfai_root/kinfai0/haorli/KaiWang_dRNA/K562_directRNA/pod5_pass/ -o pod5_pass.blow5
/scratch/kinfai_root/kinfai0/haorli/software/f5c_ReDD/f5c index pod5_pass.fastq.bgz --slow5 pod5_pass.blow5
minimap2 -ax splice -uf -k14 --secondary=no -t6 /scratch/kinfai_root/kinfai0/haorli/ReDD_rna004/reference/genome.fa.bgz pod5_pass.fastq.bgz | samtools view -h -F 2308 - | samtools sort -@ 6 - -o pod5_pass.bam && samtools index pod5_pass.bam


/scratch/kinfai_root/kinfai0/haorli/software/f5c_ReDD/f5c eventalign --rna --min-mapq 0 --pore rna004 -r pod5_pass.fastq.bgz -b pod5_pass.bam -g /scratch/kinfai_root/kinfai0/haorli/ReDD_rna004/reference/genome.fa.bgz --slow5 pod5_pass.blow5 \
-o window_size_9_hg38 \
--redd-candidate /scratch/kinfai_root/kinfai0/haorli/ReDD_rna004/reditools2_candidates/K562.candidate_sites.tab \
--redd --redd-window-size 9 --scale-events -B 20M -K 30720