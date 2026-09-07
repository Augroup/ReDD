# ReDD: **R**NA **e**diting detection by **D**irect RNA sequencing and **D**eep learning
* ReDD: **R**NA **e**diting detection by **D**irect RNA sequencing and **D**eep learning
* [Argo-ReDD](https://redd-portal.me): a **cloud-based platform** to run ReDD online without the need to prepare a computing environment.
* This repository contains code and tutorials to run ReDD.
* **RNA004 data** (pod5 + Dorado) is supported through a separate workflow, see [RNA004 support](#rna004-support); pre-built Docker images: [`tidesun/redd-rna004`](https://hub.docker.com/r/tidesun/redd-rna004).
* The code to reproduce figures and results in the manuscript is deposited in `reproduce_scripts` folder.
# Getting started
**Besides running ReDD locally, we also provides [Argo-ReDD](https://redd-portal.me): a **cloud-based platform** to run ReDD online without the need to prepare a computing environment.**

# RNA002 support
## Download and Installation
The following comand and pipeline has been tested in the following Linux systems:
* Red Hat Enterprise Linux Server release 7.9 (Maipo)
*  Ubuntu 20.04.3 LTS

**Download**: download ReDD code through github. 
```
git clone https://github.com/Augroup/ReDD && cd ReDD
```

**Installation**: install ReDD, you should install mamba or conda first. 
To install mamba, if you already have conda 4.14.0 or upper installed you can skip this step.
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
source $HOME/.bashrc
```
Then create conda environment ReDD from file environment.yaml, tested under conda 4.14.0
```
conda env create -n ReDD --file environment.yaml
```
Before running ReDD activate the ReDD environment first
```
conda activate ReDD
```

## Usage

ReDD pipeline is supported by Snakemake, refer to  [snakemake](https://snakemake.readthedocs.io/en/stable/).   for more information of Snakemake. (todo). The whole pipeline can be found in the snakefile under the ReDD-main directory.
We used script generate_script.py to generate the config.yaml file and shell file, run.pbs, for users to run ReDD in a Slurm cluster or bash. You need to set parameter --pipeline_mode to either 'cluster' or 'bash'.
ReDD pepeline provides two ways for the mapping of reads: mapping to genome and transcriptome. You need to choose one at the begining. Mapping to the transcriptome has one additional step where the genome wide annotaions will be added to the site-level and molecule-level results, which requires a reference genome and reference annotation as well.

* **Use genome as reference:**
The results will be generated in the **{output_path}/outputs** folder. 
The molecule-level result will be generated in file **{output_name}_genome.prediction.txt** 
The site-level result before any filtering will be generated in file **{output_name}_genome.site.bed**
The site-level result after filtering by predicted editing ratios, read coverage or other filters will be generated in file **{output_name}_genome.flt.tab**
The results for visualization will be generated in the **{OUTPUT_PATH}/outputs/precomputed_visualization** folder. You can uploaded the results to ReDD cloud-based server for visualization and anayses.
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
                        {'cluster' or 'bash'}. In 'cluster' mode, ReDD will be
                        run through submitting to a Slurm cluster. In 'bash'
                        mode, ReDD will run in bash command. [default:cluster]
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
The results for visualization will be generated in the **{OUTPUT_PATH}/outputs/precomputed_visualization** folder. You can uploaded the results to ReDD cloud-based server for visualization and anayses.
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
                        {'cluster' or 'bash'}. In 'cluster' mode, ReDD will be
                        run through submitting to a Slurm cluster. In 'bash'
                        mode, ReDD will run in bash command. [default:cluster]
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
mkdir ~/ReDD_data && cd ~/ReDD_data
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
cd ~/ReDD
python generate_script.py transcriptome \
--pipeline_mode 'bash' \
--output_name H1-DE_sample_trans \
--output_path ~/ReDD-results/ \
--input_fastq ~/ReDD_data/fastq_pass \
--input_fast5 ~/ReDD_data/fast5_pass \
--input_summary ~/ReDD_data/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome ~/ReDD_data/genome.fa \
--ref_transcriptome ~/ReDD_data/Stem_cell_talon.flt.bam_flt.gtf.fa \
--device 'CPU' \
--ref_annotation ~/ReDD_data/Stem_cell_talon.flt.bam_flt.gpd
```
```
cd ~/ReDD
python generate_script.py genome \
--pipeline_mode 'bash' \
--output_name H1-DE_sample_genome \
--output_path ~/ReDD-results/ \
--input_fastq ~/ReDD_data/fastq_pass \
--input_fast5 ~/ReDD_data/fast5_pass \
--input_summary ~/ReDD_data/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome ~/ReDD_data/genome.fa \
--device 'CPU' \
--ref_alu ~/ReDD_data/Hg38_Alu.merge.bed \
--ref_snp ~/ReDD_data/hg38_snp151.bed \
--ref_REDIportal ~/ReDD_data/REDIportal_hg38.txt 
```

# RNA004 support

ReDD now supports ONT **RNA004** direct RNA data (kit `SQK-RNA004`, flow cell `FLO-PRO004RA`, pod5 raw data, Dorado
basecalling). RNA004 has its own workflow ([Snakefile_RNA004](Snakefile_RNA004)), its own conda environment
([environment_rna004.yaml](environment_rna004.yaml)) and a new PyTorch model; the RNA002 pipeline described above is unchanged.

## Installation (RNA004)

```
git clone https://github.com/Augroup/ReDD && cd ReDD
conda env create -n ReDD_RNA004 --file environment_rna004.yaml   # run from the repository root
conda activate ReDD_RNA004
bash scripts/rna004/install_dorado.sh                            # Dorado 1.1.1 + rna004_130bps_sup@v5.2.0 -> software/dorado
bash scripts/rna004/download_model.sh                            # RNA004 model weights -> scripts/models/rna004/general.pt
```
The environment build compiles and installs the vendored **uncalled4 ReDD fork** in [uncalled4_ReDD/](uncalled4_ReDD/) (source of the
`--redd-out`, `--redd-candidate` and `--redd-window-size` options; see [uncalled4_ReDD/VENDORED_FROM.txt](uncalled4_ReDD/VENDORED_FROM.txt)).
PyTorch 2.6.0 is installed from PyPI (runs on CPU and CUDA GPUs). `install_dorado.sh` downloads Dorado and the basecalling model
from Oxford Nanopore into `software/dorado/` (Dorado has its own licence and is not stored in this repository); Dorado needs a GPU
for realistic data sizes.

**Model weights**: `scripts/models/rna004/general.pt` (652 MB) is not stored in git; download it before running the pipeline:
```
bash scripts/rna004/download_model.sh      # -> scripts/models/rna004/general.pt (from https://reddexamples.s3.us-east-2.amazonaws.com/RNA004/general.pt)
```

## Docker (RNA004)

Ready-to-run images are published on Docker Hub as [`tidesun/redd-rna004`](https://hub.docker.com/r/tidesun/redd-rna004):

| tag | torch build | use |
|---|---|---|
| `gpu` (= `latest`) | CUDA 12.4 wheel | NVIDIA GPU hosts (driver + [NVIDIA container toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/) needed); also runs on CPU |
| `cpu` | CPU-only wheel (smaller) | hosts without a GPU |

The image contains the pipeline code (`main` branch, at `/opt/ReDD`) and the `ReDD_RNA004` conda environment. Dorado, the basecalling
model and the ReDD model weights (~8 GB) are downloaded inside the container by `redd-fetch-assets`; run it at the start of your
command as below. Everything else in `generate_script.py rna004` then works with its default paths.

```
# your data: download assets, generate run.pbs/config.yaml, execute the pipeline
docker run --rm --gpus all -v /path/on/host:/data tidesun/redd-rna004:gpu bash -c \
  "redd-fetch-assets && \
   python generate_script.py rna004 --pipeline_mode bash --input_pod5 /data/pod5 --ref_genome /data/genome.fa \
     --output_path /data/run1 --output_name sample1 --device GPU && bash /data/run1/run.pbs"

# smoke test on the bundled test data (600 HEK293T-WT reads, from pod5 including basecalling); the test data is not in the image
git clone --depth 1 https://github.com/Augroup/ReDD.git
docker run --rm --gpus all -v $PWD/ReDD/test_data:/opt/ReDD/test_data -v /path/on/host:/data tidesun/redd-rna004:gpu bash -c \
  "redd-fetch-assets && bash test_data/rna004/run_test.sh /data/redd_test GPU"
```
Without a GPU use `tidesun/redd-rna004:cpu`, drop `--gpus all` and pass `--device CPU` (see the CPU notes below; Dorado on CPU is slow).
Already-basecalled data can be given with `--input_bam /data/dorado.bam` (Dorado BAM with move tags) instead of basecalling inside the
container. The assets are re-downloaded for every new container; to download them once and reuse them, add a volume:
`docker volume create redd-assets` and `-v redd-assets:/opt/redd-assets` on every `docker run` (`redd-fetch-assets` skips files that
are already present). `DORADO_VERSION`, `DORADO_MODEL` and `REDD_MODEL` (e.g. `-e DORADO_MODEL=rna004_130bps_hac@v5.2.0`) select other
Dorado/model versions.

With Apptainer/Singularity (e.g. on an HPC cluster; the container filesystem is read-only, so give the assets a host directory):
```
apptainer pull redd-rna004.sif docker://tidesun/redd-rna004:gpu
mkdir -p redd-assets
apptainer exec --nv -B redd-assets:/opt/redd-assets -B /path/on/host:/data redd-rna004.sif bash -c \
  "redd-fetch-assets && cd /opt/ReDD && python generate_script.py rna004 --pipeline_mode bash --input_pod5 /data/pod5 \
     --ref_genome /data/genome.fa --output_path /data/run1 --output_name sample1 --device GPU && bash /data/run1/run.pbs"
```

## Usage (RNA004)

```
python generate_script.py rna004 --help
```
Required: `--input_pod5` (pod5 file or directory), `--ref_genome`, `--output_path`, `--output_name`, `--device {CPU,GPU}`.
By default the pipeline basecalls the pod5 with Dorado (`--dorado_bin`, `--dorado_model`, `--dorado_device`, `--threads_basecall`);
to reuse existing basecalls give `--input_bam` (Dorado BAM produced with `--emit-moves`; **`--emit-fastq` output cannot be used**,
it drops the move tags) or `--input_fastq` (FASTQ made with `samtools fastq -T "mv,ts,pi,sp,ns" dorado.bam`).

Useful options: `--ref_candidate_sites` (REDItools-style table; candidate sites are labelled with their bulk editing ratio and written
to a separate feature file, both candidate and non-candidate A sites are predicted), `--num_split` (number of contig groups processed in
parallel, default 24), `--threads_extract`, `--threads_predict`, `--batch_size`, `--slurm_gpu_args` (sbatch options added to the
basecalling and prediction jobs in cluster mode when `--device GPU`, default `--gpus-per-node=1`), `--conda_env`, and the same
site-level filtering/annotation options as the RNA002 pipeline (`--ref_alu`, `--ref_snp`, `--ref_REDIportal`, `--filter_snp`,
`--filter_m6A`, `--coverage_cutoff`, `--ratio_cutoff`, ...).

`generate_script.py rna004` writes `run.pbs` and `config.yaml` into `--output_path`; run `bash run.pbs` (`--pipeline_mode bash`) or
`sbatch run.pbs` (`--pipeline_mode cluster`). Pipeline steps:

1. `dorado basecaller <model> pod5/ --emit-moves` (skipped when `--input_bam`/`--input_fastq` is given).
2. `samtools fastq -T` keeps the Dorado tags, `minimap2 -y -ax splice -uf -k14 --secondary=no` aligns to the genome, primary alignments are kept (`igv/{name}.genome.sorted.bam`).
3. The BAM is split into contig groups; each group is processed by `uncalled4 align --redd-out` (17-base windows around every reference A) into `intermediates/cache/{name}.{group}[.candidate|.noncandidate].hdf5`.
4. `scripts/rna004/predict.py` computes the per-read editing probability of every window (`intermediates/cache/{name}.{group}.prediction.raw.txt`). Predictions are per read and per window, so the split into groups does not change any value.
5. Molecule-level results are merged into **`outputs/{name}.prediction.genome.txt`**
   (`label  read_id  contig  position(1-based)  strand  probability`; label is `I:<ratio>` for candidate sites and `I` otherwise),
   site-level results into **`outputs/{name}.site.bed`** (before filtering) and **`outputs/{name}.flt.genome.tab`** (after filtering/annotation),
   and the visualization pre-computation into `outputs/precomputed_visualization/{name}/`, exactly as for the RNA002 genome pipeline.

## Reproducibility notes (RNA004)

* Feature extraction is deterministic and per read: given the same alignment, `uncalled4 align --redd-out` reproduces feature windows
  bit-for-bit, whether the signal is read from pod5 directly or from a blow5 conversion of the same pod5, and independently of how
  reads are split into groups.
* Prediction is per window (`eval()` mode: BatchNorm running statistics, no dropout), so the read set, batch size and contig
  grouping do not change any probability. Probabilities computed on a GPU and on a CPU differ slightly (median |Δp| ≈ 4e-4,
  up to a few 1e-2 for borderline windows) because cuDNN uses TF32 arithmetic by default on Ampere and newer GPUs; `predict.py`
  keeps the PyTorch defaults so GPU runs reproduce the published model outputs.
* Dorado 1.1.1 cannot open pod5 files written by pod5 library ≥ 0.3.30 (read table v6); MinKNOW output is fine, but if you subset
  or convert pod5 files yourself use `pod5 <= 0.3.28`.
* **Blackwell and newer GPUs (compute capability sm_120: RTX 50-series, RTX PRO 6000, B100/B200)**: the validated `torch==2.6.0`
  build has no kernels for them, so `predict.py` detects this and falls back to CPU with a warning (`--force_gpu` overrides). To predict on
  such a GPU, create the environment from [environment_rna004_cu128.yaml](environment_rna004_cu128.yaml) (torch 2.7.1+cu128) instead;
  probabilities may then differ from the published torch-2.6.0 outputs at the ~1e-3 level. Dorado is unaffected (it bundles its own CUDA).
* Dorado basecalling is meant for GPUs. On CPU the pipeline passes `--batchsize 8` (override with `--dorado_batchsize`): with
  Dorado's automatic batch size the `sup` model needs more than 64 GB of RAM. As a reference, 278 reads take ~12 min / 26 GB with
  `sup` and ~2 min / 14 GB with `hac` on 8 CPU threads. Dorado's CPU and GPU basecalls are not identical (sequences and move
  tables differ slightly), so probabilities are only reproducible from the same basecalls (i.e. from a given Dorado BAM) or when
  basecalling on the same hardware class.

## Test run (RNA004)

[test_data/rna004/](test_data/rna004/) contains a small set of HEK293T-WT RNA004 reads as MinKNOW pod5 (plus the same reads as a
Dorado BAM with move tags), a 360 kb reference (`chr11_sub` = hg38 chr11:600,001-960,000) and the matching candidate sites.
With the environment activated and Dorado installed:
```
bash test_data/rna004/run_test.sh [output_dir] [CPU|GPU] [pod5|bam]
```
`pod5` (default) runs the full pipeline including Dorado basecalling; `bam` starts from the provided basecalls. Outputs land in
`output_dir/outputs/`: `HEK293T-WT_chr11_sub.prediction.genome.txt` (178,594 read x A-site predictions), `.site.bed` (2,845 sites
with coverage >= 5), `.flt.genome.tab` (1,083 sites after filtering) and `precomputed_visualization/`. On 8 CPU threads the `bam` mode
takes ~40 min (uncalled4 extraction ~30 s, model inference ~37 min); on a GPU inference takes seconds. The `pod5` mode adds Dorado
`sup` basecalling of the 600 reads (~26 min and ~28 GB RAM on 8 CPU threads; seconds on a GPU) and gives 178,029 predictions /
2,823 sites / 1,068 filtered sites — slightly different from the `bam` mode because CPU and GPU basecalls of the same reads differ.

## Possible issues and solutions (to be continue)
* After run generate_script.py, a **run.pbs**, and a **config.yaml** file will be generated in {output_path}.You can refere to run.pbs and config.yaml for details of the commands and configurations.
* For Slurm submission (--pipeline_mode='cluster') , you can modify **run.pbs** according to your cluster's system configuration.
* Each time of running **run.pbs**, a **log** file and a **dag.pdf** will be generated in **{output_path}/ReDD_logs/ReDD_{output_name}.log** and **{output_path}**. You can refere to the log file for running status and dag.pdf for the rule dependencies at each step of the ReDD. The tags in the dag.pdf file represent the rule names which can be found in the Snakefile.
* **Simply rerun the pipeline by sbatch run.pbs or bash run.pbs usually solves most issues. The pipeline will continue from where it failed last time.** 
* If pipeline is stopped due to time limit for "cluster" mode, you can resubmit the job by "sbatch run.pbs" to make run for another {overall_time} or enlarge the {overall_time} parameter
* If pipeline is stopped due to memory issues, you can try resubmit the job or reduce the number of reads per task by {num_reads} parameter.
* If you want to change some parameters and rerun the pipline from where the parameter effective, you should change the parameters in the config.yaml file manualy and remove the results in the "ouputs" folder or "intermediates" folder from where the parameter effective and rerun run.pbs again. The relationship of output file can be found in the Snakemake file according to the rule dependencies presented in the dag.pdf.


### Acknowledgements
Research and development of the methods in this package were funded by
NIH grant XXXXX

### Citation
If you use any part of this code in your work, please cite our
[ReDD paper](http://).

### License
This software is released under the MIT license. For more details, please refer
[LICENSE.txt](https://github.com/Tidesun/ReDD/LICENSE.txt).

For questions, please email wangdu@missouri.edu.
