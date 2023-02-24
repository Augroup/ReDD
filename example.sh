python generate_script.py transcriptome \
--pipeline_mode 'cluster' \
--output_name H1-DE_trans_onecell \
--overall_time 24 \
--account PCON0009 \
--max_cores_resources 336 \
--output_path ./ \
--filter_snp \
--filter_m6A \
--input_fastq /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/fastq_pass \
--input_fast5 /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--ref_transcriptome /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gtf.fa \
--ref_annotation /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gpd \
--ref_alu /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/Hg38_Alu.merge.bed \
--ref_snp /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/hg38_snp151.bed \
--ref_REDIportal /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/REDIportal_hg38.txt \
--ref_candidate_sites /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/DE-H1_directRNA.candidate_sites.tab \
--device 'GPU'
#Submitted batch job 13094756


python generate_script.py genome \
--pipeline_mode 'cluster' \
--input_fastq /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/fastq_pass \
--input_fast5 /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--output_name H1-DE_genome \
--output_path /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDD-results/H1-DE_genome \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--overall_time 24 \
--account PCON0009 \
--max_cores_resources 112 \
--filter_snp \
--filter_m6A \
--ref_alu /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/Hg38_Alu.merge.bed \
--ref_snp /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/hg38_snp151.bed \
--ref_REDIportal /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/REDIportal_hg38.txt \
--ref_candidate_sites /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/DE-H1_directRNA.candidate_sites.tab \
--device 'GPU'

#13094760


python generate_script.py genome \
--pipeline_mode 'cluster' \
--input_fastq /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample2/fastq_pass \
--input_fast5 /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample2/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--output_name H1-DE_sample2 \
--output_path /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDD-results/H1-DE_sample2 \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--account PCON0009 \
--max_cores_resources 112 \
--ref_candidate_sites /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/DE-H1_directRNA.candidate_sites.tab \
--device 'GPU'



python generate_script.py genome \
--pipeline_mode 'bash' \
--input_fastq /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample2/fastq_pass \
--input_fast5 /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample2/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--output_name H1-DE_sample2 \
--output_path /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDD-results/H1-DE_sample2 \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--account PCON0009 \
--ref_candidate_sites /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/DE-H1_directRNA.candidate_sites.tab \
--device 'CPU'


python generate_script.py transcriptome \
--pipeline_mode 'bash' \
--output_name H1-DE_trans \
--overall_time 24 \
--account PCON0009 \
--output_path /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDD-results/H1-DE_sample2 \
--input_fastq /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample2/fastq_pass \
--input_fast5 /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample2/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--ref_transcriptome /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gtf.fa \
--ref_annotation /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gpd \
--ref_alu /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/Hg38_Alu.merge.bed \
--ref_snp /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/hg38_snp151.bed \
--ref_REDIportal /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/REDIportal_hg38.txt \



python generate_script.py transcriptome \
--pipeline_mode 'cluster' \
--output_name H1-DE_sample_trans \
--account PCON0009 \
--output_path /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDD-results/H1-DE_sample \
--input_fastq /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample/fastq_pass \
--input_fast5 /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--ref_transcriptome /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gtf.fa \
--device 'GPU' \
--ref_annotation /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gpd \
--ref_alu /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/Hg38_Alu.merge.bed \
--ref_snp /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/hg38_snp151.bed \
--ref_REDIportal /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/REDIportal_hg38.txt 

#Submitted batch job 13108401


python generate_script.py genome \
--pipeline_mode 'cluster' \
--output_name H1-DE_sample_genome \
--account PCON0009 \
--output_path /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDD-results/H1-DE_sample \
--input_fastq /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample/fastq_pass \
--input_fast5 /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--account PCON0009 \
--ref_candidate_sites /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/DE-H1_directRNA.candidate_sites.tab \
--device 'GPU'



python generate_script.py genome \
--pipeline_mode 'bash' \
--output_name H1-DE_sample_genome \
--account PCON0009 \
--output_path /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDD-results/H1-DE_sample_genome \
--input_fastq /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample/fastq_pass \
--input_fast5 /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--account PCON0009 \
--ref_candidate_sites /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/DE-H1_directRNA.candidate_sites.tab \
--device 'CPU'




python generate_script.py transcriptome \
--pipeline_mode 'cluster' \
--output_name H1-DE_sample_trans \
--account PCON0009 \
--output_path /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDD-results/H1-DE_sample_trans \
--input_fastq /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample/fastq_pass \
--input_fast5 /fs/ess/scratch/PCON0009/duolin/REDD-pipeline/REDDdata/H1-DE_sample/fast5_pass \
--input_summary /fs/project/PCON0009/LabData/Human/DE-H1/Transcriptome/ONT/directRNA/20210330_2124_X2_FAP47598_1c046625/sequencing_summary_FAP47598_07e34f33.txt \
--ref_genome /fs/project/PCON0009/Au-scratch2/ying/StemCell/ref/genome.fa \
--ref_transcriptome /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gtf.fa \
--device 'GPU' \
--ref_annotation /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/2-curate_gtf/Stem_cell_talon.flt.bam_flt.gpd \
--ref_alu /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/Hg38_Alu.merge.bed \
--ref_snp /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/hg38_snp151.bed \
--ref_REDIportal /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test/REDIportal_hg38.txt 
#submit 13108362 