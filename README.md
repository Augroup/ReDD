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
                                 [--ref_candidate_sites REF_CANDIDATE_SITES]

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
  --ref_candidate_sites REF_CANDIDATE_SITES
                        The path of reference candidate sites
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
                                        [--ref_candidate_sites REF_CANDIDATE_SITES]
                                        --ref_transcriptome REF_TRANSCRIPTOME
                                        --ref_annotation REF_ANNOTATION

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
  --ref_candidate_sites REF_CANDIDATE_SITES
                        The path of reference candidate sites
  --ref_transcriptome REF_TRANSCRIPTOME
                        The path of reference transcriptome
  --ref_annotation REF_ANNOTATION
                        The path of reference annotation
```
Then run the script or submit job
```
sbatch run.pbs
```
or
```
bash run.pbs
```
## Examples
```
python ~/_projects/AIediting_pipelines/generate_script.py transcriptome \
--output_name H1-DE_trans \
--overall_time 24 \
--account PCON0009 \
--max_cores_resources 336 \
--num_split 4 \
--output_path ./ \
--filter_snp True \
--filter_m6A True \
--input_fastq /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/fastq_pass \
--input_fast5 /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--ref_transcriptome /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gtf.fa \
--ref_annotation /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/gencode.v31.annotation.gpd \
--ref_alu /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/Hg38_Alu.merge.bed \
--ref_snp /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/hg38_snp151.bed \
--ref_REDIportal /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/REDIportal_hg38.txt \
--ref_candidate_sites /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/DE-H1_directRNA.candidate_sites.tab
--snakefile_path /users/PCON0009/haoranli/_projects/AIediting_pipelines/

python ~/_projects/AIediting_pipelines/generate_script.py genome \
--output_name H1-DE_genome \
--overall_time 24 \
--account PCON0009 \
--max_cores_resources 112 \
--num_split 4 \
--output_path ./ \
--filter_snp True \
--filter_m6A True \
--input_fastq /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/fastq_pass \
--input_fast5 /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--ref_alu /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/Hg38_Alu.merge.bed \
--ref_snp /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/hg38_snp151.bed \
--ref_m6A_motif /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/analysis_extract_m6A_motif.bed \
--ref_REDIportal /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/REDIportal_hg38.txt \
--ref_candidate_sites /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/DE-H1_directRNA.candidate_sites.tab
--snakefile_path /users/PCON0009/haoranli/_projects/AIediting_pipelines/
```
