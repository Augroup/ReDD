#!/bin/bash
cd /fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/test
##cdna_sites filtering pipeline
input_txt=/fs/ess/scratch/PCON0009/duolin/modification/hg38cDNAown_Eventfeature_5feature_win9/hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_mergedH_epochs35_KOnonw5withLSTM/AFG-H1_directRNA_onlycandidate.txt
input_bed=AFG-H1_directRNA_onlycandidate.cdna_model.bed
cdna_output=AFG-H1_directRNA_onlycandidate.cdna_model.cov5.anno

awk -F '\t' -v OFS="\t" '{hash_total[$3"\t"$4]++; if($6>=0.5) hash_edit[$3"\t"$4]++; hash_score[$3"\t"$4]+=$6; }END{for(site in hash_total){edit=0;if(site in hash_edit) edit=hash_edit[site]; ratio=edit/hash_total[site]; score=hash_score[site]/hash_total[site]; printf site"\t"edit"\t"hash_total[site]"\t%.3f\t%.3f\n",ratio,score; } }' $input_txt | sort -k1,1 -k2,2n > $input_bed
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$1"\t"$2]++}ARGIND==2{if($1"\t"$2 in hash_cdna){ hash_cdna_site[$1"\t"$2]=$3"\t"$4; hash_cdna_strand[$1"\t"$2]=$5; hash_site[$3"\t"$4]++}}ARGIND==3{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==4{hash_gene[$2]=$1}ARGIND==5{if($1"\t"$2 in hash_site) hash_snp[$1"\t"$2]++}ARGIND==6{if($1"\t"$2 in hash_site) hash_m6a[$1"\t"$2]++}ARGIND==7{site=hash_cdna_site[$1"\t"$2]; db="not_in_REDIportal"; snp="non_snp"; m6a="non_m6A_motif"; if(site in hash_db) db="REDIportal"; if(site in hash_snp) snp="snp"; if(site in hash_m6a) m6a="m6A_motif"; print $0,site,hash_gene[$1],hash_cdna_strand[$1"\t"$2],db,snp,m6a; }' $input_bed gencode.v31.annotation.cdna2genome.tab REDIportal_hg38.txt gencode.v31.annotation.gpd hg38_snp151.bed analysis_extract_m6A_motif.bed $input_bed > $cdna_output-1.tab
awk -F '\t' -v OFS="\t" '{print $7,$8,$8+1}' $cdna_output-1.tab | sort -k1,1 -k2,2n | uniq > $cdna_output-1.bed
bedtools intersect -a $cdna_output-1.bed -b Hg38_Alu.merge.bed -wo > $cdna_output-2.bed
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_alu[$1"\t"$2]++}ARGIND==2{site=$7"\t"$8; label="non_Alu"; gene=$9; if(site in hash_alu) label="Alu"; if($4>=5) print $0,label; }' $cdna_output-2.bed $cdna_output-1.tab > $cdna_output.tab
rm $cdna_output-1.tab
rm $cdna_output-1.bed
rm $cdna_output-2.bed
##flt_snp
awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($12=="non_snp") print $0 }' $cdna_output.tab > $cdna_output.flt_snp.tab
##flt_m6A
awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($13=="non_m6A_motif") print $0 }' $cdna_output.tab > $cdna_output.flt_m6A.tab

##cdna2genome
input_txt=/fs/ess/scratch/PCON0009/duolin/modification/hg38cDNAown_Eventfeature_5feature_win9/hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_mergedH_epochs35_KOnonw5withLSTM/AFG-H1_directRNA_onlycandidate.txt
cdna2genome_output=AFG-H1_directRNA_onlycandidate.cdna_model.cdna2genome

awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$3"\t"$4]++}ARGIND==2{if($1"\t"$2 in hash_cdna){ hash_cdna_site[$1"\t"$2]=$3"\t"$4; hash_cdna_strand[$1"\t"$2]=$5; hash_site[$3"\t"$4]++}}ARGIND==3{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==4{hash_gene[$2]=$1}ARGIND==5{site=hash_cdna_site[$3"\t"$4]; db="not_in_REDIportal"; if(site in hash_db) db="REDIportal"; print $3,$4,$2,$6,site,hash_gene[$3],hash_cdna_strand[$3"\t"$4],db; }' $input_txt gencode.v31.annotation.cdna2genome.tab REDIportal_hg38.txt gencode.v31.annotation.gpd $input_txt > $cdna2genome_output.txt
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_site[$7"\t"$8]=$11"\t"$12"\t"$13"\t"$14}ARGIND==2{hash_total[$5"\t"$6]++; if($4>=0.5) hash_edit[$5"\t"$6]++; hash_score[$5"\t"$6]+=$4; }END{for(site in hash_total){if(site in hash_site) {edit=0;if(site in hash_edit) edit=hash_edit[site]; ratio=edit/hash_total[site]; score=hash_score[site]/hash_total[site]; print site,edit,hash_total[site],ratio,score,hash_site[site]}}}' $cdna_output.tab $cdna2genome_output.txt | sort -k1,1 -k2,2n > $cdna2genome_output.tab
rm $cdna2genome_output.txt
##flt_snp
awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($8=="non_snp") print $0 }' $cdna2genome_output.tab > $cdna2genome_output.flt_snp.tab
##flt_m6A
awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($9=="non_m6A_motif") print $0 }' $cdna2genome_output.tab > $cdna2genome_output.flt_m6A.tab

##genome_sites filtering pipeline
genome_input_txt=/fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/Duolin_hg38_predict_2/AFG-H1_directRNA_onlycandidate.txt
genome_input=AFG-H1_directRNA_onlycandidate.genome_model.bed
genome_output=AFG-H1_directRNA_onlycandidate.genome_model.cov5.anno

awk -F '\t' -v OFS="\t" '{hash_total[$3"\t"$4]++; if($6>=0.5) hash_edit[$3"\t"$4]++; hash_score[$3"\t"$4]+=$6; }END{for(site in hash_total){edit=0;if(site in hash_edit) edit=hash_edit[site]; ratio=edit/hash_total[site]; score=hash_score[site]/hash_total[site]; printf site"\t"edit"\t"hash_total[site]"\t%.3f\t%.3f\n",ratio,score; } }' $genome_input_txt | sort -k1,1 -k2,2n > $genome_input
awk -F '\t' -v OFS="\t" '{print $1,$2,$2+1}' $genome_input | sort -k1,1 -k2,2n | uniq > $genome_output-1.bed
bedtools intersect -a $genome_output-1.bed -b Hg38_Alu.merge.bed -wo > $genome_output-2.bed
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_site[$1"\t"$2]++ }ARGIND==2{if($1"\t"$2 in hash_site) hash_snp[$1"\t"$2]++}ARGIND==3{if($1"\t"$2 in hash_site) hash_m6a[$1"\t"$2]++}ARGIND==4{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==5{hash_alu[$1"\t"$2]++}ARGIND==6{site=$1"\t"$2; db="not_in_REDIportal"; snp="non_snp"; m6a="non_m6A_motif"; if(site in hash_db) db="REDIportal"; if(site in hash_snp) snp="snp"; if(site in hash_m6a) m6a="m6A_motif"; label="non_Alu"; if(site in hash_alu) label="Alu"; if($4>=5) print $0,db,snp,m6a,label; }' $genome_input hg38_snp151.bed analysis_extract_m6A_motif.bed REDIportal_hg38.txt $genome_output-2.bed $genome_input > $genome_output.tab
rm $genome_output-1.bed
rm $genome_output-2.bed
##flt_snp
awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($8=="non_snp") print $0 }' $genome_output.tab > $genome_output.flt_snp.tab
##flt_m6A
awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($9=="non_m6A_motif") print $0 }' $genome_output.tab > $genome_output.flt_m6A.tab