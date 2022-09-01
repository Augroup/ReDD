cdna_molecule_input=$1
genome_molecule_output=$2
genome_site_output=$3
cdna2genome_tab=$4
REDI_txt=$5
annotation_gpd=$6
cdna_site_tab=$7
filter_snp=$8
filter_m6A=$9

awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$3"\t"$4]++}ARGIND==2{if($1"\t"$2 in hash_cdna){ hash_cdna_site[$1"\t"$2]=$3"\t"$4; hash_cdna_strand[$1"\t"$2]=$5; hash_site[$3"\t"$4]++}}ARGIND==3{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==4{hash_gene[$2]=$1}ARGIND==5{site=hash_cdna_site[$3"\t"$4]; db="not_in_REDIportal"; if(site in hash_db) db="REDIportal"; print db,$2,site,hash_cdna_strand[$3"\t"$4],$6; }' $cdna_molecule_input $cdna2genome_tab $REDI_txt $annotation_gpd $cdna_molecule_input > $genome_molecule_output
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$3"\t"$4]++}ARGIND==2{if($1"\t"$2 in hash_cdna){ hash_cdna_site[$1"\t"$2]=$3"\t"$4; hash_cdna_strand[$1"\t"$2]=$5; hash_site[$3"\t"$4]++}}ARGIND==3{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==4{hash_gene[$2]=$1}ARGIND==5{site=hash_cdna_site[$3"\t"$4]; db="not_in_REDIportal"; if(site in hash_db) db="REDIportal"; print $3,$4,$2,$6,site,hash_gene[$3],hash_cdna_strand[$3"\t"$4],db; }' $cdna_molecule_input $cdna2genome_tab $REDI_txt $annotation_gpd $cdna_molecule_input > $genome_molecule_output.temp
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_site[$7"\t"$8]=$11"\t"$12"\t"$13"\t"$14}ARGIND==2{hash_total[$5"\t"$6]++; if($4>=0.5) hash_edit[$5"\t"$6]++; hash_score[$5"\t"$6]+=$4; }END{for(site in hash_total){if(site in hash_site) {edit=0;if(site in hash_edit) edit=hash_edit[site]; ratio=edit/hash_total[site]; score=hash_score[site]/hash_total[site]; print site,edit,hash_total[site],ratio,score,hash_site[site]}}}' $cdna_site_tab $genome_molecule_output.temp | sort -k1,1 -k2,2n > $genome_site_output
# rm $genome_molecule_output.temp
if [ $filter_snp = True ];
then
    mv $genome_site_output $genome_site_output.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($8=="non_snp") print $0 }' $genome_site_output.temp > $genome_site_output
    rm $genome_site_output.temp
fi
if [ $filter_m6A = True ];
then
    mv $genome_site_output $genome_site_output.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($9=="non_m6A_motif") print $0 }' $genome_site_output.temp > $genome_site_output
    rm $genome_site_output.temp
fi