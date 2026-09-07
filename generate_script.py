import argparse
from pathlib import Path
import os
import numpy as np
#from datetime import datetime
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="REDD pipelines script generator",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")
    genome_parser = subparsers.add_parser('genome',help='Use genome as reference')
    transcriptome_parser = subparsers.add_parser('transcriptome',help='Use transcriptome as reference')
    add_rna004_parser(subparsers)
    
    g_required_group = genome_parser.add_argument_group('Required arguments about input,output and reference')
    t_required_group = transcriptome_parser.add_argument_group('Required arguments about input,output and reference')
    g_pipeline_group = genome_parser.add_argument_group('Optional arguments about other pipeline settings')
    t_pipeline_group = transcriptome_parser.add_argument_group('Optional arguments about other pipeline settings')
    g_ref_group = genome_parser.add_argument_group('Optional arguments about reference files')
    t_ref_group = transcriptome_parser.add_argument_group('Optional arguments about reference files')
    for re_group,name in zip([g_required_group,t_required_group],['genome_parser','transcriptome_parser']):
        #input
        re_group.add_argument('--input_fastq', type=str, help="The path of directory containing input fastq files (/fastq_pass/)",required=True)
        re_group.add_argument('--input_fast5', type=str, help="The path of directory containing input fast5 files  (/fast5_pass/)",required=True)
        re_group.add_argument('--input_summary', type=str, help="The path of sequencing summary file (sequencing_summary_*.txt)",required=True)
        #output
        re_group.add_argument('-o','--output_path', type=str, help="The path of output directory, multiple runs can share a [output_path], the specific run is indicated by [output_name].",required=True)
        re_group.add_argument('-n','--output_name', type=str, help="The name to output." ,required=True)
        #reference
        re_group.add_argument('--ref_genome', type=str, help="The path of reference genome",required=True)
        if name == 'transcriptome_parser':
            re_group.add_argument('--ref_transcriptome', type=str, help="The path of reference transcriptome",required=True)
            # ref_group.add_argument('--ref_cdna_to_genome', type=str, help="The path of reference cdna to genome file",required=True)
            re_group.add_argument('--ref_annotation', type=str, help="The path of reference annotation",required=True)
        
        #device
        re_group.add_argument('--device',type=str,help="Use GPU or CPU to speed up prediction select from {'CPU','GPU'}.[default:GPU]",required=True)
    
    for pipeline_group in [g_pipeline_group,t_pipeline_group]:
        # pipeline groups
        pipeline_group.add_argument('--overall_time', type=int, help="The estimated overall running time. It will only affect the runing time of one-time submittion. [default:24]",default=24,required=False)
        pipeline_group.add_argument('--pipeline_mode', type=str, help="Which mode to use to run the pipeline. Choose from {'cluster' or 'bash'}. In 'cluster' mode, REDD will be run through submitting to a Slurm cluster. In 'bash' mode, REDD will run in bash command. [default:cluster]",default='cluster',required=False)
        pipeline_group.add_argument('--account', type=str, help="Which account uses for running in cluster. Required if --pipeline_mode set to 'cluster'",default=None,required=False)
        pipeline_group.add_argument('--max_cores_resources', type=int, help="Max number of cores uses in the same time.[default:112]",default=112,required=False)
        pipeline_group.add_argument('--num_reads', type=int, help="Number of reads per task to enable parallel computing.[default:100,000]. Only valid, if --num_split is not set.",default=100000,required=False)
        pipeline_group.add_argument('--num_split', type=int, help="Number of split task to enable parallel computing.[default: None]. If not set, it will be calculated by parameter --num_reads. If set, --num_reads is invalid.",default=None,required=False)
        pipeline_group.add_argument('--model',type=str,help="Which model will be used {'general','stemcell_specific','HEK293T_specific','GM12878_specific'}[default:general]",default='general',required=False)
        pipeline_group.add_argument('--filter_snp',action='store_true',help="Once set, the SNP sites in the SNP annotation provided by --ref_snp will be removed from site-level results.",default='False',required=False)
        pipeline_group.add_argument('--filter_m6A',action='store_true',help="Whether filer out sites with m6A motifs. If set, sites that have m6A motifs will be generated from the reference genome and removed from site-level output.",required=False)
        
        pipeline_group.add_argument('--coverage_cutoff', type=int, help="read coverage threshod for site level outputs, if --ref_alu is not provided.[default:10]",default=10,required=False)
        pipeline_group.add_argument('--ratio_cutoff', type=float, help="ratio threshod for site level outputs,if --ref_alu is not provided.[default:0.1]",default=0.1,required=False)
    
    # reference group
    for ref_group,name in zip([g_ref_group,t_ref_group],['genome_parser','transcriptome_parser']):
        ref_group.add_argument('--ref_alu', type=str, help="The path of reference ALU. Refer to ./scripts/reference/Hg38_Alu.merge.bed for a format example. If provided, the site-level results will be filtered based on ALU regions.",default=None,required=False)
        ref_group.add_argument('--in_alu_ratio', type=float, help="ratio_cutoff for sites in ALU region. Only valid when --ref_alu is provided [default:0.1]",default=0.1,required=False)
        ref_group.add_argument('--out_alu_ratio', type=float, help="ratio_cutoff for sites not in ALU region.Only valid when --ref_alu is provided [default:0.3]",default=0.3,required=False)
        ref_group.add_argument('--in_alu_coverage', type=int, help="threshod for read coverage in ALU region. Only valid when --ref_alu is provided.[default:10]",default=10,required=False)
        ref_group.add_argument('--out_alu_coverage', type=int, help="threshod for read coverage not in ALU region. Only valid when --ref_alu is provided.[default:30]",default=30,required=False)
        ref_group.add_argument('--ref_snp', type=str, help="The path of reference SNP, if provided the SNP information will be added to the site-level results and filtering by SNPs is available.",default=None,required=False)
        ref_group.add_argument('--add_m6A_reference', help="If set, the m6A motif information will be calculated and added to the site-level results based on the reference genome.",action='store_true',required=False)
        ref_group.add_argument('--ref_REDIportal', type=str, help="The path of reference REDIportal,if provided the annotation in REDIportal will be added to the site-level results",required=False,default=None)
        ref_group.add_argument('--ref_candidate_sites', type=str, help="The path of reference candidate sites",required=False,default=None)
    
    args = parser.parse_args()
    return args

def prompt_error(error_str):
    print(error_str)
    exit()

def read_counts(input_fastq_folder):
    count=0
    
    for fastqfile in os.listdir(input_fastq_folder):
      
      if fastqfile.endswith(".fastq"):
        
        file=open(input_fastq_folder+"/"+fastqfile)
        for line in file:
            count+=1
    
    return int(count/4)



# ----------------------------------------------------------------------------------------------
# RNA004 (pod5 + Dorado BAM, uncalled4 feature extraction, torch model). Uses Snakefile_RNA004.
# ----------------------------------------------------------------------------------------------
RNA004_WINDOW_SIZE = 17

def add_rna004_parser(subparsers):
    p = subparsers.add_parser('rna004', help='RNA004 (SQK-RNA004 / FLO-PRO004RA) direct RNA data, genome reference')
    req = p.add_argument_group('Required arguments about input, output and reference')
    req.add_argument('--input_pod5', type=str, required=True, help="pod5 file, or directory containing pod5 files (searched recursively). slow5/blow5/fast5 are accepted as well.")
    reads = req.add_mutually_exclusive_group(required=False)
    reads.add_argument('--input_bam', type=str, help="Already basecalled reads from Dorado as (unaligned) BAM with move tags (dorado basecaller ... --emit-moves). A file or a directory of BAM files. If neither --input_bam nor --input_fastq is given, the pipeline basecalls the pod5 files with Dorado (see --dorado_*).")
    reads.add_argument('--input_fastq', type=str, help="Already basecalled reads as FASTQ whose headers carry the Dorado move tags (made with `samtools fastq -T 'mv,ts,pi,sp,ns' dorado.bam`). A file or directory. Plain FASTQ without mv/ts tags cannot be used.")
    req.add_argument('-o', '--output_path', type=str, required=True, help="The path of output directory, multiple runs can share a [output_path], the specific run is indicated by [output_name].")
    req.add_argument('-n', '--output_name', type=str, required=True, help="The name to output.")
    req.add_argument('--ref_genome', type=str, required=True, help="The path of reference genome (FASTA)")
    req.add_argument('--device', type=str, required=True, help="Use GPU or CPU for prediction, select from {'CPU','GPU'}")

    pipe = p.add_argument_group('Optional arguments about other pipeline settings')
    pipe.add_argument('--overall_time', type=int, default=24, help="Estimated overall running time (hours) of the one-time submission. [default:24]")
    pipe.add_argument('--pipeline_mode', type=str, default='cluster', help="{'cluster','bash'}. 'cluster' submits every step to a Slurm cluster, 'bash' runs locally. [default:cluster]")
    pipe.add_argument('--account', type=str, default=None, help="Slurm account. Required if --pipeline_mode is 'cluster'")
    pipe.add_argument('--max_cores_resources', type=int, default=112, help="Max number of cores used at the same time. [default:112]")
    pipe.add_argument('--num_split', type=int, default=24, help="Number of reference-contig groups processed in parallel (feature extraction and prediction). [default:24, or the number of contigs if smaller]")
    pipe.add_argument('--threads_extract', type=int, default=8, help="Threads per uncalled4 feature-extraction job. [default:8]")
    pipe.add_argument('--threads_predict', type=int, default=4, help="CPU threads per prediction job. [default:4]")
    pipe.add_argument('--batch_size', type=int, default=4096, help="Prediction batch size. [default:4096]")
    pipe.add_argument('--bam_chunksize', type=int, default=100, help="uncalled4 --bam-chunksize. [default:100]")
    pipe.add_argument('--model', type=str, default='general', help="Model name under scripts/models/rna004/ (without .pt), or a path to a checkpoint. [default:general]")
    pipe.add_argument('--flowcell', type=str, default='FLO-PRO004RA', help="uncalled4 --flowcell [default:FLO-PRO004RA]")
    pipe.add_argument('--kit', type=str, default='SQK-RNA004', help="uncalled4 --kit [default:SQK-RNA004]")
    pipe.add_argument('--dorado_bin', type=str, default=None, help="Dorado executable used to basecall the pod5 files. [default: <ReDD>/software/dorado/bin/dorado, installed by scripts/rna004/install_dorado.sh]")
    pipe.add_argument('--dorado_model', type=str, default='rna004_130bps_sup@v5.2.0', help="Dorado basecalling model: a directory, or a model name looked up in <ReDD>/software/dorado/models. [default:rna004_130bps_sup@v5.2.0]")
    pipe.add_argument('--dorado_mod_models', type=str, default=None, help="Optional comma-separated modified-base models passed to dorado --modified-bases-models (does not change the canonical basecalls or move tags).")
    pipe.add_argument('--dorado_device', type=str, default=None, help="dorado -x device. [default: 'cuda:all' if --device GPU, else 'cpu']")
    pipe.add_argument('--dorado_batchsize', type=int, default=None, help="dorado --batchsize. [default: Dorado auto-selection on GPU; 8 on CPU, where the auto-selected batch needs >64 GB RAM for the sup model]")
    pipe.add_argument('--threads_basecall', type=int, default=8, help="CPU threads for the Dorado basecalling job. [default:8]")
    pipe.add_argument('--conda_env', type=str, default='ReDD_RNA004', help="Conda environment name or path activated in run.pbs. [default:ReDD_RNA004]")
    pipe.add_argument('--slurm_gpu_args', type=str, default='--gpus-per-node=1', help="Extra sbatch arguments added to prediction jobs when --device GPU in cluster mode, e.g. '--gpus-per-node=1 --partition=gpu'. [default:--gpus-per-node=1]")
    pipe.add_argument('--filter_snp', action='store_true', help="Remove SNP sites (from --ref_snp) from site-level results.")
    pipe.add_argument('--filter_m6A', action='store_true', help="Remove sites with m6A motifs (derived from the reference genome) from site-level results.")
    pipe.add_argument('--coverage_cutoff', type=int, default=10, help="Read coverage threshold for site-level outputs, if --ref_alu is not provided. [default:10]")
    pipe.add_argument('--ratio_cutoff', type=float, default=0.1, help="Ratio threshold for site-level outputs, if --ref_alu is not provided. [default:0.1]")

    ref = p.add_argument_group('Optional arguments about reference files')
    ref.add_argument('--ref_candidate_sites', type=str, default=None, help="Candidate editing sites (REDItools-style tab: contig, 1-based pos, change, ratio, coverage). If given, uncalled4 labels candidate sites with their ratio and writes candidate/noncandidate feature files separately; both are predicted.")
    ref.add_argument('--ref_annotation', type=str, default=None, help="Reference annotation in GPD format, used for visualization pre-computation.")
    ref.add_argument('--ref_alu', type=str, default=None, help="Reference ALU bed; if provided, site-level results are filtered by ALU regions.")
    ref.add_argument('--in_alu_ratio', type=float, default=0.1)
    ref.add_argument('--out_alu_ratio', type=float, default=0.3)
    ref.add_argument('--in_alu_coverage', type=int, default=10)
    ref.add_argument('--out_alu_coverage', type=int, default=30)
    ref.add_argument('--ref_snp', type=str, default=None, help="Reference SNP bed.")
    ref.add_argument('--add_m6A_reference', action='store_true', help="Add m6A motif information to site-level results.")
    ref.add_argument('--ref_REDIportal', type=str, default=None, help="REDIportal annotation.")
    return p


def contig_lengths(ref_genome):
    fai = ref_genome + '.fai'
    if not os.path.exists(fai):
        try:
            import pysam
            pysam.faidx(ref_genome)
        except Exception as e:
            prompt_error(f'{fai} not found and could not be created ({e}). Run: samtools faidx {ref_genome}')
    lengths = []
    with open(fai) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 2:
                lengths.append((fields[0], int(fields[1])))
    if not lengths:
        prompt_error(f'No contigs found in {fai}')
    return lengths


def make_region_groups(lengths, num_split):
    """Greedy balanced partition of contigs into num_split groups by length."""
    num_split = max(1, min(num_split, len(lengths)))
    groups = [[] for _ in range(num_split)]
    totals = [0] * num_split
    for contig, length in sorted(lengths, key=lambda x: -x[1]):
        i = totals.index(min(totals))
        groups[i].append(contig)
        totals[i] += length
    width = len(str(num_split))
    return {f'g{str(i + 1).zfill(width)}': g for i, g in enumerate(groups) if g}


def yaml_list(items):
    return '[' + ', '.join(f"'{x}'" for x in items) + ']'


def main_rna004(args):
    output_folder = os.path.abspath(args.output_path)
    output_name = args.output_name
    snakefile_path = os.path.dirname(os.path.realpath(__file__))
    if args.pipeline_mode == 'cluster' and args.account is None:
        prompt_error('No account(--account) is given but you have set --pipeline_mode to cluster')
    if args.filter_snp and args.ref_snp is None:
        prompt_error('No reference snp file(--ref_snp) is given but you have set --filter_snp')
    if args.device not in ('CPU', 'GPU'):
        prompt_error("--device must be 'CPU' or 'GPU'")

    if os.path.exists(args.model):
        model_path = os.path.abspath(args.model)
    else:
        model_path = os.path.join(snakefile_path, 'scripts', 'models', 'rna004', args.model + '.pt')
    if not os.path.exists(model_path):
        prompt_error(f'Model checkpoint not found: {model_path}')

    input_pod5 = os.path.abspath(args.input_pod5)
    dorado_bin = dorado_model = ''
    if args.input_bam is not None:
        input_type, input_reads = 'bam', os.path.abspath(args.input_bam)
    elif args.input_fastq is not None:
        input_type, input_reads = 'fastq', os.path.abspath(args.input_fastq)
    else:
        input_type, input_reads = 'pod5', input_pod5
        dorado_bin = args.dorado_bin or os.path.join(snakefile_path, 'software', 'dorado', 'bin', 'dorado')
        if not os.path.exists(dorado_bin):
            prompt_error(f'Dorado not found at {dorado_bin}. Install it with `bash scripts/rna004/install_dorado.sh` or pass --dorado_bin / --input_bam.')
        dorado_bin = os.path.abspath(dorado_bin)
        if os.path.isdir(args.dorado_model):
            dorado_model = os.path.abspath(args.dorado_model)
        else:
            dorado_model = os.path.join(snakefile_path, 'software', 'dorado', 'models', args.dorado_model)
            if not os.path.isdir(dorado_model):
                prompt_error(f'Dorado model not found: {dorado_model}. Run `bash scripts/rna004/install_dorado.sh` or pass a model directory to --dorado_model.')
    dorado_device = args.dorado_device or ('cuda:all' if args.device == 'GPU' else 'cpu')
    dorado_mod_models = '' if args.dorado_mod_models is None else args.dorado_mod_models
    if args.dorado_batchsize is not None:
        dorado_batchsize = args.dorado_batchsize
    else:
        dorado_batchsize = 0 if dorado_device.startswith('cuda') else 8   # 0 = let Dorado choose
    ref_genome_file = os.path.abspath(args.ref_genome)
    for path in (input_pod5, input_reads, ref_genome_file):
        if not os.path.exists(path):
            prompt_error(f'Input not found: {path}')

    Path(output_folder).mkdir(exist_ok=True, parents=True)
    Path(os.path.join(output_folder, 'REDD_logs')).mkdir(exist_ok=True, parents=True)
    empty_file = f'{output_folder}/empty'
    with open(empty_file, 'w'):
        pass

    def opt(path):
        return empty_file if path is None else os.path.abspath(path)

    region_groups = make_region_groups(contig_lengths(ref_genome_file), args.num_split)
    print(f'{len(region_groups)} contig groups will be processed in parallel.')
    region_yaml = '\n'.join(f"  {name}: {yaml_list(contigs)}" for name, contigs in region_groups.items())

    config = f"""
pore: 'rna004'
reference: 'genome'
sample: '{output_name}'
input_type: '{input_type}'
input_reads: '{input_reads}'
input_pod5: '{input_pod5}'
dorado_bin: '{dorado_bin}'
dorado_model: '{dorado_model}'
dorado_mod_models: '{dorado_mod_models}'
dorado_device: '{dorado_device}'
dorado_batchsize: {dorado_batchsize}
threads_basecall: {args.threads_basecall}
ref_genome_file: '{ref_genome_file}'
ref_candidate_sites_file: '{opt(args.ref_candidate_sites)}'
has_candidate: {args.ref_candidate_sites is not None}
ref_annotation_file: '{opt(args.ref_annotation)}'
ref_alu_file: '{opt(args.ref_alu)}'
ref_snp_file: '{opt(args.ref_snp)}'
ref_REDIportal_file: '{opt(args.ref_REDIportal)}'
empty_file: '{empty_file}'
window_size: {RNA004_WINDOW_SIZE}
flowcell: '{args.flowcell}'
kit: '{args.kit}'
bam_chunksize: {args.bam_chunksize}
threads_extract: {args.threads_extract}
threads_predict: {args.threads_predict}
device: '{args.device}'
batch_size: {args.batch_size}
model_path: '{model_path}'
slurm_gpu_args: '{args.slurm_gpu_args}'
coverage_cutoff: '{args.coverage_cutoff}'
ratio_cutoff: '{args.ratio_cutoff}'
in_alu_coverage: '{args.in_alu_coverage}'
out_alu_coverage: '{args.out_alu_coverage}'
in_alu_ratio: '{args.in_alu_ratio}'
out_alu_ratio: '{args.out_alu_ratio}'
add_m6A_reference: '{args.add_m6A_reference}'
filter_snp: '{args.filter_snp}'
filter_m6A: '{args.filter_m6A}'
region_groups:
{region_yaml}
"""
    if args.device == 'GPU':
        sbatch_line = "#SBATCH --nodes=1 --ntasks-per-node=1 --gpus-per-node=1"
    else:
        sbatch_line = "#SBATCH --nodes=1 --ntasks-per-node=1"
    sbatch_account = f"#SBATCH --account={args.account}" if args.account else ""
    if args.pipeline_mode == 'cluster':
        cluster = (f'--cluster "sbatch -A {args.account} -t {{resources.runtime}} --mem {{resources.mem_mb}} -N 1 -c {{threads}} '
                   f'{{resources.slurm_extra}} -o {output_folder}/REDD_logs/slurm-%j.out" --jobs {max(4, len(region_groups))} ')
    else:
        cluster = f'--jobs {max(4, len(region_groups))} '
    script = f"""#!/bin/bash
#SBATCH --time={args.overall_time}:00:00
{sbatch_line}
{sbatch_account}

source activate {args.conda_env}
cd {output_folder}
mkdir -p intermediates/reference intermediates/fastq intermediates/reads intermediates/mapped_reads intermediates/cache outputs igv REDD_logs
# link reference
ln -sf {ref_genome_file} intermediates/reference/genome.fa
if [ -f {ref_genome_file}.fai ]; then ln -sf {ref_genome_file}.fai intermediates/reference/genome.fa.fai; fi

# copy pipeline code (models are referenced in place)
mkdir -p scripts
cp {snakefile_path}/scripts/*.py {snakefile_path}/scripts/*.sh {snakefile_path}/scripts/*.pl scripts/
cp -r {snakefile_path}/scripts/rna004 scripts/
cp {snakefile_path}/Snakefile_RNA004 Snakefile
snakemake --unlock
snakemake -p --rulegraph | dot -Tpdf > dag.pdf

echo "Please refer to log file in {output_folder}/REDD_logs/REDD_{output_name}.log for more information."
snakemake -c {args.max_cores_resources} -p {cluster}--latency-wait 60 --rerun-incomplete --keep-going > REDD_logs/REDD_{output_name}.log 2>&1
"""
    with open(f'{output_folder}/run.pbs', 'w') as f:
        f.write(script)
    with open(f'{output_folder}/config.yaml', 'w') as f:
        f.write(config)
    print(f'run.pbs and config.yaml have been generated in {output_folder}')
    print("To run the ReDD RNA004 pipeline:\n")
    print("cd " + output_folder)
    if args.pipeline_mode == 'cluster':
        print("sbatch run.pbs")
    else:
        print("bash run.pbs")



def main():
    args = parse_arguments()
    if args.subparser_name == 'rna004':
        main_rna004(args)
        return
    output_folder = args.output_path
    output_folder = os.path.abspath(output_folder)
    output_name = args.output_name
    #timestr=datetime.now().strftime('%Y-%m-%d-%H-%M')
    input_fastq_folder = args.input_fastq
    input_fast5_folder = args.input_fast5
    input_summary_file = args.input_summary
    reference_choice = args.subparser_name
    if args.subparser_name == 'transcriptome':
        ref_transcriptome_file = args.ref_transcriptome
        # ref_cdna_to_genome_file = args.ref_cdna_to_genome
        ref_annotation_file = args.ref_annotation
    elif args.subparser_name == 'genome':
        ref_transcriptome_file = None
        # ref_cdna_to_genome_file = None
        ref_annotation_file = None
    else:
        prompt_error(f'Invalid subcommand {args.subparser_name} given!')
    
    ref_genome_file = args.ref_genome
    ref_alu_file = args.ref_alu
    ref_snp_file = args.ref_snp
    # ref_m6A_motif_file = args.ref_m6A_motif
    ref_REDIportal_file = args.ref_REDIportal
    ref_candidate_sites_file = args.ref_candidate_sites

    overall_time = args.overall_time
    account = args.account
    max_cores_resources = args.max_cores_resources
    model = args.model
    coverage_cutoff=args.coverage_cutoff
    ratio_cutoff=args.ratio_cutoff
    in_alu_coverage=args.in_alu_coverage
    out_alu_coverage=args.out_alu_coverage
    in_alu_ratio=args.in_alu_ratio
    out_alu_ratio=args.out_alu_ratio
    add_m6A_reference=args.add_m6A_reference
    filter_snp = args.filter_snp
    filter_m6A = args.filter_m6A
    device = args.device
    pipeline_mode = args.pipeline_mode
    snakefile_path = os.path.dirname(os.path.realpath(__file__)) #os.path.abspath("./")
    Path(output_folder).mkdir(exist_ok=True,parents=True)
    Path(os.path.abspath(output_folder+"/REDD_logs")).mkdir(exist_ok=True,parents=True)
    # error handling
    
    if filter_snp == 'True' and ref_snp_file is None:
        prompt_error('No reference snp file(--ref_snp) is given but you have set --filter_snp to True')
    
    if pipeline_mode == 'cluster' and account is None:
        prompt_error('No account(--account) is given but you have set --pipeline_mode to cluster')
    with open(f'{output_folder}/empty', 'w') as fp:
        pass
    if reference_choice == 'genome':
        #os.system(f"attrib +h {output_folder}/empty")
        ref_transcriptome_file = f'{output_folder}/empty'
        # ref_cdna_to_genome_file = f'{output_folder}/empty'
        ref_annotation_file = f'{output_folder}/empty'
        ref_dump_position = 'disk'
    else:
        if ref_transcriptome_file is None:
            prompt_error('No reference transcriptome(--ref_transcriptome)  is given but you have set --reference_choice to transcriptome')
        # if ref_cdna_to_genome_file is None:
        #     prompt_error('No reference cdna_to_genome_file(--ref_cdna_to_genome)  is given but you have set --reference_choice to transcriptome')
        if ref_annotation_file is None:
            prompt_error('No reference annotation(--ref_annotation) is given but you have set --reference_choice to transcriptome')
        ref_dump_position = 'memory'
    
    ref_alu_file = f'{output_folder}/empty' if ref_alu_file is None else ref_alu_file
    ref_snp_file = f'{output_folder}/empty' if ref_snp_file is None else ref_snp_file
    ref_REDIportal_file = f'{output_folder}/empty' if ref_REDIportal_file is None else ref_REDIportal_file
    ref_candidate_sites_file = f'{output_folder}/empty' if ref_candidate_sites_file is None else ref_candidate_sites_file
    num_split = args.num_split
    num_reads=args.num_reads
    if num_split is None:
        totalreads=read_counts(input_fastq_folder)
        print("Because --num_split is not set, it will be calculated based on total number of reads in your fastq folder:"+str(input_fastq_folder))
        print("total number of reads="+str(totalreads))
        num_split=int(np.ceil(totalreads/num_reads))
        print("Task will be splitted into "+str(num_split)+" tasks for parallel computting.")
    
    # ref_m6A_motif_file = f'{output_folder}/empty' if ref_m6A_motif_file is None else ref_m6A_motif_file
    
    config = f'''
reference: '{reference_choice}' # genome or transcriptome
center: 'A'
nt: 4
featurenum: 5
buffersize: 1000
labeltype: 'I'
ref_dump_position: "{ref_dump_position}"
device: '{device}'
batch_size: 10000
model: '{model}'
coverage_cutoff: '{coverage_cutoff}'
ratio_cutoff: '{ratio_cutoff}'
in_alu_coverage: '{in_alu_coverage}'
out_alu_coverage: '{out_alu_coverage}'
in_alu_ratio: '{in_alu_ratio}'
out_alu_ratio: '{out_alu_ratio}'
add_m6A_reference: '{add_m6A_reference}'
filter_snp: '{filter_snp}'
filter_m6A: '{filter_m6A}'
sample: '{output_name}'
input_fastq_folder: '{input_fastq_folder}'
input_fast5_folder: '{input_fast5_folder}'
input_summary_file: '{input_summary_file}'
ref_transcriptome_file: '{ref_transcriptome_file}'
ref_genome_file: '{ref_genome_file}'
ref_annotation_file: '{ref_annotation_file}'
ref_alu_file: '{ref_alu_file}'
ref_snp_file: '{ref_snp_file}'
ref_REDIportal_file: '{ref_REDIportal_file}'
ref_candidate_sites_file: '{ref_candidate_sites_file}'
'''
    if device=='GPU':
       SBATCHline="#SBATCH --nodes=1 --ntasks-per-node=1 --gpus-per-node=1"
    else:
       SBATCHline="#SBATCH --nodes=1 --ntasks-per-node=1"
    script = f'''#!/bin/bash
#SBATCH --time={overall_time}:00:00
'''+f'''{SBATCHline}
'''+f'''
#SBATCH --account={account}

source activate REDD
mkdir {output_folder}/intermediates
mkdir {output_folder}/intermediates/fastq/
mkdir {output_folder}/intermediates/fast5/
mkdir {output_folder}/intermediates/reference/
# link reference
ln -s {ref_transcriptome_file} {output_folder}/intermediates/reference/transcriptome.fa
ln -s {ref_genome_file} {output_folder}/intermediates/reference/genome.fa

# link script
cp {snakefile_path}/scripts -r scripts
cp {snakefile_path}/Snakefile Snakefile
snakemake --unlock
snakemake -p outputs/precomputed_visualization/{output_name} --rulegraph | dot -Tpdf > dag.pdf

echo "Please refer to log file in {output_folder}/REDD_logs/REDD_{output_name}.log for more information."
'''

    if pipeline_mode == 'cluster':
       scriptsnakemake = f'''
snakemake -c {max_cores_resources} -p outputs/precomputed_visualization/{output_name} --cluster "sbatch -A {account} -t {{resources.runtime}} --mem {{resources.mem_mb}} -N 1 -c {{threads}}" --jobs {num_split} --set-scatter split={num_split} --latency-wait 60 --rerun-incomplete > REDD_logs/REDD_{output_name}.log 2>&1
    '''
    elif pipeline_mode =='bash':
         scriptsnakemake = f'''
snakemake -c {max_cores_resources} -p outputs/precomputed_visualization/{output_name}  --jobs {num_split} --set-scatter split={num_split} --latency-wait 60 --rerun-incomplete > REDD_logs/REDD_{output_name}.log 2>&1
    '''
    
    with open(f'{output_folder}/run.pbs','w') as f:
        f.write(script)
        f.write("\n")
        f.write(scriptsnakemake)
    with open(f'{output_folder}/config.yaml','w') as f:
        f.write(config)
    
    print(f'run_{output_name}.pbs and your current configuaration config.yaml have been generated in {output_folder}')
    print("To run REDD pipeline:\n")
    print("cd "+output_folder)
    if pipeline_mode == 'cluster':
          print(f'sbatch run.pbs\n\nsince you set the pipeline_mode to \'cluster\'.')
    elif pipeline_mode == "bash":
            print(f'bash run.pbs\n\nsince you set the pipeline_mode to \'bash\'.')
    

main()
    
    
    
    
    
    
