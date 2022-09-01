import argparse
from pathlib import Path
import os
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="REDD pipelines script generator",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")
    genome_parser = subparsers.add_parser('genome',help='Use genome as reference')
    transcriptome_parser = subparsers.add_parser('transcriptome',help='Use transcriptome as reference')
    
    g_output_group = genome_parser.add_argument_group('Named arguments about outputs')
    t_output_group = transcriptome_parser.add_argument_group('Named arguments about outputs')
    g_input_group = genome_parser.add_argument_group('Named arguments about inputs')
    t_input_group = transcriptome_parser.add_argument_group('Named arguments about inputs')
    g_pipeline_group = genome_parser.add_argument_group('Named arguments about other pipeline settings')
    t_pipeline_group = transcriptome_parser.add_argument_group('Named arguments about other pipeline settings')
    g_ref_group = genome_parser.add_argument_group('Named arguments about reference files')
    t_ref_group = transcriptome_parser.add_argument_group('Named arguments about reference files')
    for output_group in [g_output_group,t_output_group]:
        # output group
        output_group.add_argument('-n','--output_name', type=str, help="The name of output results and results will be at OUTPUT_PATH/outputs/OUTPUT_NAME.prediction.txt",required=True)
        output_group.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    for input_group in [g_input_group,t_input_group]:
        # input group
        input_group.add_argument('--input_fastq', type=str, help="The path of directory containing input fastq files (/fastq_pass/)",required=True)
        input_group.add_argument('--input_fast5', type=str, help="The path of directory containing input fast5 files  (/fast5_pass/)",required=True)
        input_group.add_argument('--input_summary', type=str, help="The path of sequencing summary file (sequencing_summary_*.txt)",required=True)
    for pipeline_group in [g_pipeline_group,t_pipeline_group]:
        # pipeline groups
        pipeline_group.add_argument('--overall_time', type=int, help="The estimated overall running time.[default:24]",default=24,required=False)
        pipeline_group.add_argument('--pipeline_mode', type=str, help="Which mode to use to run the pipeline.[cluster,bash][default:cluster]",default='cluster',required=False)
        pipeline_group.add_argument('--account', type=str, help="Which account uses for running in cluster. Required if --pipeline_mode set to cluster",default=None,required=False)
        pipeline_group.add_argument('--max_cores_resources', type=int, help="Max number of cores uses in the same time.[default:112]",default=112,required=False)
        pipeline_group.add_argument('--num_split', type=int, help="Number of file split to enable parallel computing.[default:4]",default=4,required=False)
        pipeline_group.add_argument('--model',type=str,help="Which model will be used [general][default:general]",default='general',required=False)
        pipeline_group.add_argument('--filter_snp',type=str,help="Whether filer out snp[default:False]",default='False',required=False)
        pipeline_group.add_argument('--filter_m6A',type=str,help="Whether filer out m6A[default:False]",default='False',required=False)
        pipeline_group.add_argument('--device',type=str,help="Use GPU or CPU to speed up prediction [CPU,GPU].[default:CPU]",default='CPU',required=False)
        pipeline_group.add_argument('--snakefile_path',type=str,help="The path of directory containing snakefile.[default:./]",default='./',required=False)
    
    # reference group
    for ref_group,name in zip([g_ref_group,t_ref_group],['genome_parser','transcriptome_parser']):
        ref_group.add_argument('--ref_genome', type=str, help="The path of reference genome",required=True)
        ref_group.add_argument('--ref_alu', type=str, help="The path of reference alu",required=True)
        ref_group.add_argument('--ref_snp', type=str, help="The path of reference snp",default=None,required=False)
        ref_group.add_argument('--ref_m6A_motif', type=str, help="The path of reference m6A motif",default=None,required=False)
        ref_group.add_argument('--ref_REDIportal', type=str, help="The path of reference REDI portal file",required=True)
        ref_group.add_argument('--ref_candidate_sites', type=str, help="The path of reference candidate sites",required=False,default=None)
        if name == 'transcriptome_parser':
            ref_group.add_argument('--ref_transcriptome', type=str, help="The path of reference transcriptome",required=True)
            # ref_group.add_argument('--ref_cdna_to_genome', type=str, help="The path of reference cdna to genome file",required=True)
            ref_group.add_argument('--ref_annotation', type=str, help="The path of reference annotation",required=True)

    args = parser.parse_args()
    return args
def prompt_error(error_str):
    print(error_str)
    exit()


def main():
    args = parse_arguments()
    output_folder = args.output_path
    output_folder = os.path.abspath(output_folder)
    output_name = args.output_name
    
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
    num_split = args.num_split
    model = args.model
    filter_snp = args.filter_snp
    filter_m6A = args.filter_m6A
    device = args.device
    pipeline_mode = args.pipeline_mode
    snakefile_path = args.snakefile_path
    Path(output_folder).mkdir(exist_ok=True,parents=True)
    # error handling
    
    if filter_snp == 'True' and ref_snp_file is None:
        prompt_error('No reference snp file(--ref_snp) is given but you have set --filter_snp to True')
    # if filter_m6A == 'True' and ref_m6A_motif_file is None:
    #     prompt_error('No reference m6A file(--ref_m6A_motif) is given but you have set --filter_m6A to True')
    
    if pipeline_mode == 'cluster' and account is None:
        prompt_error('No account(--account) is given but you have set --pipeline_mode to cluster')
    if reference_choice == 'genome':
        with open(f'{output_folder}/empty', 'w') as fp:
            pass
        os.system(f"attrib +h {output_folder}/empty")
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
    ref_snp_file = f'{output_folder}/empty' if ref_alu_file is None else ref_alu_file
    ref_candidate_sites_file = f'{output_folder}/empty' if ref_candidate_sites_file is None else ref_candidate_sites_file
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
coverage_cutoff: 5
ratio_cutoff: 0
filter_snp: '{filter_snp}'
filter_m6A: '{filter_m6A}'
'''
    
    script = f'''#!/bin/bash
#SBATCH --time={overall_time}:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --account={account}

source activate AIediting

mkdir {output_folder}/intermediates
mkdir {output_folder}/intermediates/fastq/
mkdir {output_folder}/intermediates/fast5/
mkdir {output_folder}/intermediates/summary/
mkdir {output_folder}/intermediates/reference/
# link data
ln -s {input_fastq_folder} {output_folder}/intermediates/fastq/{output_name}
ln -s {input_fast5_folder} {output_folder}/intermediates/fast5/{output_name}
ln -s {input_summary_file} {output_folder}/intermediates/summary/{output_name}_summary.txt
# link reference
ln -s {ref_transcriptome_file} {output_folder}/intermediates/reference/transcriptome.fa
ln -s {ref_genome_file} {output_folder}/intermediates/reference/genome.fa
ln -s {ref_annotation_file} {output_folder}/intermediates/reference/gencode.v31.annotation.gpd
ln -s {ref_alu_file} {output_folder}/intermediates/reference/Hg38_Alu.merge.bed
ln -s {ref_snp_file} {output_folder}/intermediates/reference/hg38_snp151.bed
ln -s {ref_REDIportal_file} {output_folder}/intermediates/reference/REDIportal_hg38.txt
ln -s {ref_candidate_sites_file} {output_folder}/intermediates/reference/{output_name}.candidate_sites.tab

# link script
ln -s {snakefile_path}/scripts scripts
ln -s {snakefile_path}/Snakefile Snakefile

snakemake --unlock
snakemake -c {max_cores_resources} -p outputs/precomputed_visualization/{output_name} --cluster "sbatch -A {account} -t {{resources.runtime}} -N 1 -c {{threads}}" --jobs {num_split} --set-scatter split={num_split} --latency-wait 60 --rerun-incomplete
'''
    with open(f'{output_folder}/run.pbs','w') as f:
        f.write(script)
    with open(f'{output_folder}/config.yaml','w') as f:
        f.write(config)
main()
    
    
    
    
    
    