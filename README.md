# Installation
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
mamba env create -n AIediting --file environment.yaml
```
# Usage
```
source activate AIediting
mkdir intermediates
cd intermediates
mkdir fastq/
mkdir fast5/
mkdir summary/
mkdir reference/
ln -s /fs/project/PCON0009/LabData/Human/H1-hESC/Transcriptome/ONT/directRNA/GT20-02070_HSC-1_RNA_F_R1/fastq_pass fastq/H1-hESC
ln -s /fs/project/PCON0009/LabData/Human/H1-hESC/Transcriptome/ONT/directRNA/GT20-02070_HSC-1_RNA_F_R1/fast5_pass fast5/H1-hESC
ln -s /fs/project/PCON0009/LabData/Human/H1-hESC/Transcriptome/ONT/directRNA/GT20-02070_HSC-1_RNA_F_R1/sequencing_summary_FAL42026_033d8d92.txt summary/H1-hESC_summary.txt
ln -s /fs/project/PCON0009/Au-scratch2/haoran/reference/genecode/GRCh38.primary_assembly.genome.fa reference/genome.fa
cd ../

snakemake -c 336 -p outputs/H1-hESC_fast5_list.txt --cluster "sbatch -A PCON0009 -t {resources.runtime} -N 1 -c {threads}" --jobs 12 --set-scatter split=12
```
snakemake -c (total_number_of_CPUs[12 * 28]) -p (outputs) --cluster "sbatch -A PCON0009 -t {resources.runtime} -N 1 -c {threads}" --jobs (number_of_parallel_jobs) --set-scatter split=(number_of_samples_to_split)