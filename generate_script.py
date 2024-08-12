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



def main():
    args = parse_arguments()
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
ln -s {input_fastq_folder} {output_folder}/intermediates/fastq/{output_name}
ln -s {input_fast5_folder} {output_folder}/intermediates/fast5/{output_name}
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
    
    
    
    
    
    
