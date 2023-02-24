# REDD 
* **RNA editing using deep-learning and direct RNA sequencing revealing A-to-I editing profile during stem cell differentiation.** 

### This repository contains code and tutorials to run REDD pipeline locally
## Getting started
#### Download and Installation
The following comand and pipeline has been tested in the following Linux systems:
* Red Hat Enterprise Linux Server release 7.9 (Maipo)
*  Ubuntu 20.04.3 LTS

**Download**: download REDD code through github. 
```
git clone https://github.com/duolinwang/REDD && cd REDD
```

**Installation**: install REDD, you should install mamba or conda first. 
To install mamba, if you already have conda 4.14.0 or upper installed you can skip this step.
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
source $HOME/.bashrc
```
Then create conda environment REDD from file environment.yaml, tested under conda 4.14.0
```
conda env create -n REDD --file environment.yaml
```
Before running REDD activate the REDD environment first
```
conda activate REDD
```

# Usage

REDD pipeline is supported by Snakemake, refer to  [snakemake](https://snakemake.readthedocs.io/en/stable/).   for more information of Snakemake. (todo). The whole pipeline can be found in the snakefile under the REDD-main directory.
We used script generate_script.py to generate the config.yaml file and shell file, run.pbs, for users to run REDD in a Slurm cluster or bash. You need to set parameter --pipeline_mode to either 'cluster' or 'bash'.
REDD pepeline provides two ways for the mapping of reads: mapping to genome and transcriptome. You need to choose one at the begining. Mapping to the transcriptome has one additional step where the genome wide annotaions will be added to the site-level and molecule-level results, which requires a reference genome and reference annotation as well.

* **Use genome as reference:**
The results will be generated in the **{output_path}/outputs** folder. 
The molecule-level result will be generated in file **{output_name}_genome.prediction.txt** 
The site-level result before any filtering will be generated in file **{output_name}_genome.site.bed**
The site-level result after filtering by predicted editing ratios, read coverage or other filters will be generated in file **{output_name}_genome.flt.tab**
The results for visualization will be generated in the **{OUTPUT_PATH}/outputs/precomputed_visualization** folder. You can uploaded the results to REDD cloud-based server for visualization and anayses.
All the other intermediate results will be generated in the **intermediates** folder.

```
python generate_script.py genome --help
usage: generate_script.py transcriptome [-h] --input_fastq INPUT_FASTQ
                                        --input_fast5 INPUT_FAST5
                                        --input_summary INPUT_SUMMARY -o
                                        OUTPUT_PATH -n OUTPUT_NAME
                                        --ref_genome REF_GENOME
                                        --ref_transcriptome REF_TRANSCRIPTOME
                                        --ref_annotation REF_ANNOTATION
                                        --device DEVICE
                                        [--overall_time OVERALL_TIME]
                                        [--pipeline_mode PIPELINE_MODE]
                                        [--account ACCOUNT]
                                        [--max_cores_resources MAX_CORES_RESOURCES]
                                        [--num_reads NUM_READS]
                                        [--num_split NUM_SPLIT]
                                        [--model MODEL] [--filter_snp]
                                        [--filter_m6A]
                                        [--coverage_cutoff COVERAGE_CUTOFF]
                                        [--ratio_cutoff RATIO_CUTOFF]
                                        [--ref_alu REF_ALU]
                                        [--in_alu_ratio IN_ALU_RATIO]
                                        [--out_alu_ratio OUT_ALU_RATIO]
                                        [--in_alu_coverage IN_ALU_COVERAGE]
                                        [--out_alu_coverage OUT_ALU_COVERAGE]
                                        [--ref_snp REF_SNP]
                                        [--add_m6A_reference]
                                        [--ref_REDIportal REF_REDIPORTAL]
                                        [--ref_candidate_sites REF_CANDIDATE_SITES]

optional arguments:
  -h, --help            show this help message and exit

Required arguments about input,output and reference:
  --input_fastq INPUT_FASTQ
                        The path of directory containing input fastq files
                        (/fastq_pass/)
  --input_fast5 INPUT_FAST5
                        The path of directory containing input fast5 files
                        (/fast5_pass/)
  --input_summary INPUT_SUMMARY
                        The path of sequencing summary file
                        (sequencing_summary_*.txt)
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory, multiple runs can share
                        a [output_path], the specific run is indicated by
                        [output_name].
  -n OUTPUT_NAME, --output_name OUTPUT_NAME
                        The name to output.
  --ref_genome REF_GENOME
                        The path of reference genome
  --ref_transcriptome REF_TRANSCRIPTOME
                        The path of reference transcriptome
  --ref_annotation REF_ANNOTATION
                        The path of reference annotation
  --device DEVICE       Use GPU or CPU to speed up prediction select from
                        {'CPU','GPU'}.[default:GPU]

Optional arguments about other pipeline settings:
  --overall_time OVERALL_TIME
                        The estimated overall running time. It will only
                        affect the runing time of one-time submittion.
                        [default:24]
  --pipeline_mode PIPELINE_MODE
                        Which mode to use to run the pipeline. Choose from
                        {'cluster' or 'bash'}. In 'cluster' mode, REDD will be
                        run through submitting to a Slurm cluster. In 'bash'
                        mode, REDD will run in bash command. [default:cluster]
  --account ACCOUNT     Which account uses for running in cluster. Required if
                        --pipeline_mode set to cluster
  --max_cores_resources MAX_CORES_RESOURCES
                        Max number of cores uses in the same
                        time.[default:112]
  --num_reads NUM_READS
                        Number of reads per task to enable parallel
                        computing.[default:300,000]. Only valid, if
                        --num_split is not set.
  --num_split NUM_SPLIT
                        Number of split task to enable parallel
                        computing.[default: None]. If not set, it will be
                        calculated by parameter --num_reads. If set,
                        --num_reads is invalid.
  --model MODEL         Which model will be used {'general','stemcell_specific
                        ','HEK293T_specific','GM12878_specific'}[default:gener
                        al]
  --filter_snp          Once set, the SNP sites in the SNP annotation provided
                        by --ref_snp will be removed from site-level results.
  --filter_m6A          Whether filer out sites with m6A motifs. If set, sites
                        that have m6A motifs will be generated from the
                        reference genome and removed from site-level output.
  --coverage_cutoff COVERAGE_CUTOFF
                        read coverage threshod for site level outputs, if
                        --ref_alu is not provided.[default:10]
  --ratio_cutoff RATIO_CUTOFF
                        ratio threshod for site level outputs,if --ref_alu is
                        not provided.[default:0.1]

Optional arguments about reference files:
  --ref_alu REF_ALU     The path of reference ALU. Refer to
                        ./scripts/reference/Hg38_Alu.merge.bed for a format
                        example. If provided, the site-level results will be
                        filtered based on ALU regions.
  --in_alu_ratio IN_ALU_RATIO
                        ratio_cutoff for sites in ALU region. Only valid when
                        --ref_alu is provided [default:0.1]
  --out_alu_ratio OUT_ALU_RATIO
                        ratio_cutoff for sites not in ALU region.Only valid
                        when --ref_alu is provided [default:0.3]
  --in_alu_coverage IN_ALU_COVERAGE
                        threshod for read coverage in ALU region. Only valid
                        when --ref_alu is provided.[default:10]
  --out_alu_coverage OUT_ALU_COVERAGE
                        threshod for read coverage not in ALU region. Only
                        valid when --ref_alu is provided.[default:30]
  --ref_snp REF_SNP     The path of reference SNP, if provided the SNP
                        information will be added to the site-level results
                        and filtering by SNPs is available.
  --add_m6A_reference   If set, the m6A motif information will be calculated
                        and added to the site-level results based on the
                        reference genome.
  --ref_REDIportal REF_REDIPORTAL
                        The path of reference REDIportal,if provided the
                        annotation in REDIportal will be added to the site-
                        level results
  --ref_candidate_sites REF_CANDIDATE_SITES
                        The path of reference candidate sites
```
* **Use transcriptome as reference:**
The results will be generated in the **{output_path}/outputs** folder. 
The molecule-level result will be generated in in file **{output_name}.prediction.transcriptome.txt** 
The site-level result before any filtering will be generated in file **{output_name}.site.bed**
The site-level result after filtering by predicted editing ratios, read coverage or other filters will be generated in file **{output_name}.flt.transcriptome.tab**
The results for visualization will be generated in the **{OUTPUT_PATH}/outputs/precomputed_visualization** folder. You can uploaded the results to REDD cloud-based server for visualization and anayses.
All the other intermediate results will be generated in the **intermediates** folder.

```
python generate_script.py transcriptome --help
usage: generate_script.py transcriptome [-h] --input_fastq INPUT_FASTQ
                                        --input_fast5 INPUT_FAST5
                                        --input_summary INPUT_SUMMARY -o
                                        OUTPUT_PATH -n OUTPUT_NAME
                                        --ref_genome REF_GENOME
                                        --ref_transcriptome REF_TRANSCRIPTOME
                                        --ref_annotation REF_ANNOTATION
                                        --device DEVICE
                                        [--overall_time OVERALL_TIME]
                                        [--pipeline_mode PIPELINE_MODE]
                                        [--account ACCOUNT]
                                        [--max_cores_resources MAX_CORES_RESOURCES]
                                        [--num_reads NUM_READS]
                                        [--num_split NUM_SPLIT]
                                        [--model MODEL] [--filter_snp]
                                        [--filter_m6A]
                                        [--coverage_cutoff COVERAGE_CUTOFF]
                                        [--ratio_cutoff RATIO_CUTOFF]
                                        [--ref_alu REF_ALU]
                                        [--in_alu_ratio IN_ALU_RATIO]
                                        [--out_alu_ratio OUT_ALU_RATIO]
                                        [--in_alu_coverage IN_ALU_COVERAGE]
                                        [--out_alu_coverage OUT_ALU_COVERAGE]
                                        [--ref_snp REF_SNP]
                                        [--add_m6A_reference]
                                        [--ref_REDIportal REF_REDIPORTAL]
                                        [--ref_candidate_sites REF_CANDIDATE_SITES]

optional arguments:
  -h, --help            show this help message and exit

Required arguments about input,output and reference:
  --input_fastq INPUT_FASTQ
                        The path of directory containing input fastq files
                        (/fastq_pass/)
  --input_fast5 INPUT_FAST5
                        The path of directory containing input fast5 files
                        (/fast5_pass/)
  --input_summary INPUT_SUMMARY
                        The path of sequencing summary file
                        (sequencing_summary_*.txt)
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory, multiple runs can share
                        a [output_path], the specific run is indicated by
                        [output_name].
  -n OUTPUT_NAME, --output_name OUTPUT_NAME
                        The name to output.
  --ref_genome REF_GENOME
                        The path of reference genome
  --ref_transcriptome REF_TRANSCRIPTOME
                        The path of reference transcriptome
  --ref_annotation REF_ANNOTATION
                        The path of reference annotation
  --device DEVICE       Use GPU or CPU to speed up prediction select from
                        {'CPU','GPU'}.[default:GPU]

Optional arguments about other pipeline settings:
  --overall_time OVERALL_TIME
                        The estimated overall running time. It will only
                        affect the runing time of one-time submittion.
                        [default:24]
  --pipeline_mode PIPELINE_MODE
                        Which mode to use to run the pipeline. Choose from
                        {'cluster' or 'bash'}. In 'cluster' mode, REDD will be
                        run through submitting to a Slurm cluster. In 'bash'
                        mode, REDD will run in bash command. [default:cluster]
  --account ACCOUNT     Which account uses for running in cluster. Required if
                        --pipeline_mode set to cluster
  --max_cores_resources MAX_CORES_RESOURCES
                        Max number of cores uses in the same
                        time.[default:112]
  --num_reads NUM_READS
                        Number of reads per task to enable parallel
                        computing.[default:300,000]. Only valid, if
                        --num_split is not set.
  --num_split NUM_SPLIT
                        Number of split task to enable parallel
                        computing.[default: None]. If not set, it will be
                        calculated by parameter --num_reads. If set,
                        --num_reads is invalid.
  --model MODEL         Which model will be used {'general','stemcell_specific
                        ','HEK293T_specific','GM12878_specific'}[default:gener
                        al]
  --filter_snp          Once set, the SNP sites in the SNP annotation provided
                        by --ref_snp will be removed from site-level results.
  --filter_m6A          Whether filer out sites with m6A motifs. If set, sites
                        that have m6A motifs will be generated from the
                        reference genome and removed from site-level output.
  --coverage_cutoff COVERAGE_CUTOFF
                        read coverage threshod for site level outputs, if
                        --ref_alu is not provided.[default:10]
  --ratio_cutoff RATIO_CUTOFF
                        ratio threshod for site level outputs,if --ref_alu is
                        not provided.[default:0.1]

Optional arguments about reference files:
  --ref_alu REF_ALU     The path of reference ALU. Refer to
                        ./scripts/reference/Hg38_Alu.merge.bed for a format
                        example. If provided, the site-level results will be
                        filtered based on ALU regions.
  --in_alu_ratio IN_ALU_RATIO
                        ratio_cutoff for sites in ALU region. Only valid when
                        --ref_alu is provided [default:0.1]
  --out_alu_ratio OUT_ALU_RATIO
                        ratio_cutoff for sites not in ALU region.Only valid
                        when --ref_alu is provided [default:0.3]
  --in_alu_coverage IN_ALU_COVERAGE
                        threshod for read coverage in ALU region. Only valid
                        when --ref_alu is provided.[default:10]
  --out_alu_coverage OUT_ALU_COVERAGE
                        threshod for read coverage not in ALU region. Only
                        valid when --ref_alu is provided.[default:30]
  --ref_snp REF_SNP     The path of reference SNP, if provided the SNP
                        information will be added to the site-level results
                        and filtering by SNPs is available.
  --add_m6A_reference   If set, the m6A motif information will be calculated
                        and added to the site-level results based on the
                        reference genome.
  --ref_REDIportal REF_REDIPORTAL
                        The path of reference REDIportal,if provided the
                        annotation in REDIportal will be added to the site-
                        level results
  --ref_candidate_sites REF_CANDIDATE_SITES
                        The path of reference candidate sites
```
run.pbs and config.yaml will be generated in your output directory {OUTPUT_PATH}
Then run the script or submit job in your output directory
```
sbatch run.pbs
```
or
```
bash run.pbs
```
## Examples

**Download the test data and references from our server to a working folder**
```
mkdir ~/REDD_data && cd ~/REDD_data
```
**Download required test data**: 3959 reads randomly selected from stem cell H1-DE after basecalling
```
#Fast5 files:
wget https://reddexamples.s3.us-east-2.amazonaws.com/fast5_pass.zip
#Fastq files:
wget https://reddexamples.s3.us-east-2.amazonaws.com/fastq_pass.zip
#Sequencing summary:
wget https://reddexamples.s3.us-east-2.amazonaws.com/sequencing_summary_FAP47598_07e34f33.txt
#Required reference:
#If choose to map to genome, download the refence genome sequences.
wget https://reddexamples.s3.us-east-2.amazonaws.com/genome.fa
#If choose to map to transcriptome, download our denovo reference transcriptome for stem cell and its corresponding annotation. 
wget https://reddexamples.s3.us-east-2.amazonaws.com/Stem_cell_talon.flt.bam_flt.gtf.fa
wget https://reddexamples.s3.us-east-2.amazonaws.com/Stem_cell_talon.flt.bam_flt.gpd
```
**Other reference transcriptome and annotation can be downloaded by:**
```
#Gencode release 31 (GRCh38.p12)
wget https://reddexamples.s3.us-east-2.amazonaws.com/gencode.v31.annotation.gpd
#Our denovo reference transcriptome and annotation for GM12878 cells:
wget https://reddexamples.s3.us-east-2.amazonaws.com/GM12878_talon.flt.bam_flt.gtf.fa
wget https://reddexamples.s3.us-east-2.amazonaws.com/GM12878_talon.flt.bam_flt.gpd
#Our denovo reference transcriptome and annotation for HEK293T cells:
wget https://reddexamples.s3.us-east-2.amazonaws.com/HEK293T_talon.flt.bam_flt.gtf.fa
wget https://reddexamples.s3.us-east-2.amazonaws.com/HEK293T_talon.flt.bam_flt.gpd
#To download some optional reference files
#Reference ALUs for human hg38. If provided, the site-level results will be filtered based on ALU regions:
wget https://reddexamples.s3.us-east-2.amazonaws.com/Hg38_Alu.merge.bed
#Reference SNPs for human hg38. If provided the SNP information will be added to the site-level results and filtering by SNPs is available.
wget https://reddexamples.s3.us-east-2.amazonaws.com/hg38_snp151.bed
#Reference REDIportal for human hg38,if provided the annotation in REDIportal will be added to the site-level results.
wget https://reddexamples.s3.us-east-2.amazonaws.com/REDIportal_hg38.txt
#Candidate sites detected by bulk data
wget https://reddexamples.s3.us-east-2.amazonaws.com/DE-H1_directRNA.candidate_sites.tab
```
### Example tasks:
```
cd ~/REDD
python generate_script.py transcriptome \
--pipeline_mode 'bash' \
--output_name H1-DE_sample_trans \
--output_path ~/REDD-results/ \
--input_fastq ~/REDD_data/fastq_pass \
--input_fast5 ~/REDD_data/fast5_pass \
--input_summary ~/REDD_data/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome ~/REDD_data/genome.fa \
--ref_transcriptome ~/REDD_data/Stem_cell_talon.flt.bam_flt.gtf.fa \
--device 'CPU' \
--ref_annotation ~/REDD_data/Stem_cell_talon.flt.bam_flt.gpd
```
```
cd ~/REDD
python generate_script.py genome \
--pipeline_mode 'bash' \
--output_name H1-DE_sample_genome \
--output_path ~/REDD-results/ \
--input_fastq ~/REDD_data/fastq_pass \
--input_fast5 ~/REDD_data/fast5_pass \
--input_summary ~/REDD_data/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome ~/REDD_data/genome.fa \
--device 'CPU' \
--ref_alu ~/REDD_data/Hg38_Alu.merge.bed \
--ref_snp ~/REDD_data/hg38_snp151.bed \
--ref_REDIportal ~/REDD_data/REDIportal_hg38.txt 
```

## Possible issues and solutions (to be continue)
* After run generate_script.py, a **run.pbs**, and a **config.yaml** file will be generated in {output_path}.You can refere to run.pbs and config.yaml for details of the commands and configurations.
* For Slurm submission (--pipeline_mode='cluster') , you can modify **run.pbs** according to your cluster's system configuration.
* Each time of running **run.pbs**, a **log** file and a **dag.pdf** will be generated in **{output_path}/REDD_logs/REDD_{output_name}.log** and **{output_path}**. You can refere to the log file for running status and dag.pdf for the rule dependencies at each step of the REDD. The tags in the dag.pdf file represent the rule names which can be found in the Snakefile.
* **Simply rerun the pipeline by sbatch run.pbs or bash run.pbs usually solves most issues. The pipeline will continue from where it failed last time.** 
* If pipeline is stopped due to time limit for "cluster" mode, you can resubmit the job by "sbatch run.pbs" to make run for another {overall_time} or enlarge the {overall_time} parameter
* If pipeline is stopped due to memory issues, you can try resubmit the job or reduce the number of reads per task by {num_reads} parameter.
* If you want to change some parameters and rerun the pipline from where the parameter effective, you should change the parameters in the config.yaml file manualy and remove the results in the "ouputs" folder or "intermediates" folder from where the parameter effective and rerun run.pbs again. The relationship of output file can be found in the Snakemake file according to the rule dependencies presented in the dag.pdf.


### Acknowledgements
Research and development of the methods in this package were funded by
NIH grant XXXXX

### Citation
If you use any part of this code in your work, please cite our
[REDD paper](http://).

### License
This software is released under the MIT license. For more details, please refer
[LICENSE.txt](https://github.com/Tidesun/REDD/LICENSE.txt).

For questions, please email wangdu@missouri.edu.
