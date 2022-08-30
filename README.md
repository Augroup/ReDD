# Installation
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
mamba env create -n AIediting --file environment.yaml
```
# Usage
Use genome as reference
```
usage: generate_script.py genome [-h] -n OUTPUT_NAME -o OUTPUT_PATH
                                 --input_fastq INPUT_FASTQ --input_fast5
                                 INPUT_FAST5 --input_summary INPUT_SUMMARY
                                 [--overall_time OVERALL_TIME]
                                 [--pipeline_mode PIPELINE_MODE]
                                 [--account ACCOUNT]
                                 [--max_cores_resources MAX_CORES_RESOURCES]
                                 [--num_split NUM_SPLIT] [--model MODEL]
                                 [--filter_snp FILTER_SNP]
                                 [--filter_m6A FILTER_M6A] [--device DEVICE]
                                 [--snakefile_path SNAKEFILE_PATH]
                                 --ref_genome REF_GENOME --ref_alu REF_ALU
                                 [--ref_snp REF_SNP]
                                 [--ref_m6A_motif REF_M6A_MOTIF]
                                 --ref_REDIportal REF_REDIPORTAL

optional arguments:
  -h, --help            show this help message and exit

Named arguments about outputs:
  -n OUTPUT_NAME, --output_name OUTPUT_NAME
                        The name of output results and results will be at
                        OUTPUT_PATH/outputs/OUTPUT_NAME.prediction.txt
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory

Named arguments about inputs:
  --input_fastq INPUT_FASTQ
                        The path of directory containing input fastq files
                        (/fastq_pass/)
  --input_fast5 INPUT_FAST5
                        The path of directory containing input fast5 files
                        (/fast5_pass/)
  --input_summary INPUT_SUMMARY
                        The path of sequencing summary file
                        (sequencing_summary_*.txt)

Named arguments about other pipeline settings:
  --overall_time OVERALL_TIME
                        The estimated overall running time.[default:24]
  --pipeline_mode PIPELINE_MODE
                        Which mode to use to run the
                        pipeline.[cluster,bash][default:cluster]
  --account ACCOUNT     Which account uses for running in cluster. Required if
                        --pipeline_mode set to cluster
  --max_cores_resources MAX_CORES_RESOURCES
                        Max number of cores uses in the same
                        time.[default:112]
  --num_split NUM_SPLIT
                        Number of file split to enable parallel
                        computing.[default:4]
  --model MODEL         Which model will be used [general][default:general]
  --filter_snp FILTER_SNP
                        Whether filer out snp[default:False]
  --filter_m6A FILTER_M6A
                        Whether filer out m6A[default:False]
  --device DEVICE       Use GPU or CPU to speed up prediction
                        [CPU,GPU].[default:CPU]
  --snakefile_path SNAKEFILE_PATH
                        The path of directory containing
                        snakefile.[default:./]

Named arguments about reference files:
  --ref_genome REF_GENOME
                        The path of reference genome
  --ref_alu REF_ALU     The path of reference alu
  --ref_snp REF_SNP     The path of reference snp
  --ref_m6A_motif REF_M6A_MOTIF
                        The path of reference m6A motif
  --ref_REDIportal REF_REDIPORTAL
                        The path of reference REDI portal file
```
Use transcriptome as reference
```
usage: generate_script.py transcriptome [-h] -n OUTPUT_NAME -o OUTPUT_PATH
                                        --input_fastq INPUT_FASTQ
                                        --input_fast5 INPUT_FAST5
                                        --input_summary INPUT_SUMMARY
                                        [--overall_time OVERALL_TIME]
                                        [--pipeline_mode PIPELINE_MODE]
                                        [--account ACCOUNT]
                                        [--max_cores_resources MAX_CORES_RESOURCES]
                                        [--num_split NUM_SPLIT]
                                        [--model MODEL]
                                        [--filter_snp FILTER_SNP]
                                        [--filter_m6A FILTER_M6A]
                                        [--device DEVICE]
                                        [--snakefile_path SNAKEFILE_PATH]
                                        --ref_genome REF_GENOME --ref_alu
                                        REF_ALU [--ref_snp REF_SNP]
                                        [--ref_m6A_motif REF_M6A_MOTIF]
                                        --ref_REDIportal REF_REDIPORTAL
                                        --ref_transcriptome REF_TRANSCRIPTOME
                                        --ref_cdna_to_genome
                                        REF_CDNA_TO_GENOME --ref_annotation
                                        REF_ANNOTATION

optional arguments:
  -h, --help            show this help message and exit

Named arguments about outputs:
  -n OUTPUT_NAME, --output_name OUTPUT_NAME
                        The name of output results and results will be at
                        OUTPUT_PATH/outputs/OUTPUT_NAME.prediction.txt
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory

Named arguments about inputs:
  --input_fastq INPUT_FASTQ
                        The path of directory containing input fastq files
                        (/fastq_pass/)
  --input_fast5 INPUT_FAST5
                        The path of directory containing input fast5 files
                        (/fast5_pass/)
  --input_summary INPUT_SUMMARY
                        The path of sequencing summary file
                        (sequencing_summary_*.txt)

Named arguments about other pipeline settings:
  --overall_time OVERALL_TIME
                        The estimated overall running time.[default:24]
  --pipeline_mode PIPELINE_MODE
                        Which mode to use to run the
                        pipeline.[cluster,bash][default:cluster]
  --account ACCOUNT     Which account uses for running in cluster. Required if
                        --pipeline_mode set to cluster
  --max_cores_resources MAX_CORES_RESOURCES
                        Max number of cores uses in the same
                        time.[default:112]
  --num_split NUM_SPLIT
                        Number of file split to enable parallel
                        computing.[default:4]
  --model MODEL         Which model will be used [general][default:general]
  --filter_snp FILTER_SNP
                        Whether filer out snp[default:False]
  --filter_m6A FILTER_M6A
                        Whether filer out m6A[default:False]
  --device DEVICE       Use GPU or CPU to speed up prediction
                        [CPU,GPU].[default:CPU]
  --snakefile_path SNAKEFILE_PATH
                        The path of directory containing
                        snakefile.[default:./]

Named arguments about reference files:
  --ref_genome REF_GENOME
                        The path of reference genome
  --ref_alu REF_ALU     The path of reference alu
  --ref_snp REF_SNP     The path of reference snp
  --ref_m6A_motif REF_M6A_MOTIF
                        The path of reference m6A motif
  --ref_REDIportal REF_REDIPORTAL
                        The path of reference REDI portal file
  --ref_transcriptome REF_TRANSCRIPTOME
                        The path of reference transcriptome
  --ref_cdna_to_genome REF_CDNA_TO_GENOME
                        The path of reference cdna to genome file
  --ref_annotation REF_ANNOTATION
                        The path of reference annotation
```
