site_bed=$1
alu_bed=$2
snp_bed=$3
m6A_motif_bed=$4
REDI_txt=$5
temp_bed_1=$6
temp_bed_2=$7
site_tab=$8
awk -F '\t' -v OFS="\t" '{print $1,$2,$2+1}' $site_bed | sort -k1,1 -k2,2n | uniq > $temp_bed_1
bedtools intersect -a $temp_bed_1 -b $alu_bed -wo > $temp_bed_2
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_site[$1"\t"$2]++ }ARGIND==2{if($1"\t"$2 in hash_site) hash_snp[$1"\t"$2]++}ARGIND==3{if($1"\t"$2 in hash_site) hash_m6a[$1"\t"$2]++}ARGIND==4{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==5{hash_alu[$1"\t"$2]++}ARGIND==6{site=$1"\t"$2; if(!(site in hash_snp) && !(site in hash_m6a)){ db="not_in_REDIportal"; if(site in hash_db) db="REDIportal"; label="non_Alu"; if(site in hash_alu) label="Alu"; print $0,db,label; } }' $site_bed $snp_bed $m6A_motif_bed $REDI_txt $temp_bed_2 $site_bed > $site_tab