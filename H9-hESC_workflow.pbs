#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account=PCON0009
module load python/3.6-conda5.2
source activate AIediting

mkdir intermediates
cd intermediates
mkdir fastq/
mkdir fast5/
mkdir summary/
mkdir reference/
ln -s /fs/project/PCON0009/LabData/Human/H9-hESC/Transcriptome/ONT/directRNA/20210108_1718_X5_FAO97046_21a85c27/fastq_pass fastq/H9-hESC
ln -s /fs/project/PCON0009/LabData/Human/H9-hESC/Transcriptome/ONT/directRNA/20210108_1718_X5_FAO97046_21a85c27/fast5_pass fast5/H9-hESC
ln -s /fs/project/PCON0009/LabData/Human/H9-hESC/Transcriptome/ONT/directRNA/20210108_1718_X5_FAO97046_21a85c27/sequencing_summary_FAO97046_264259c2.txt summary/H9-hESC_summary.txt
cd ../


snakemake -c 9 -p outputs/precomputed_visualization/H9-hESC --cluster "sbatch -A PCON0009 -t {resources.runtime} -N 1 -c {threads}" --jobs 3 --set-scatter split=3 --latency-wait 60
