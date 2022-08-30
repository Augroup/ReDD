site_bed=$1
alu_bed=$2
snp_bed=$3
m6A_motif_bed=$4
REDI_txt=$5
temp_bed_1=$6
temp_bed_2=$7
site_tab=$8
filter_snp=$9
filter_m6A=${10}

awk -F '\t' -v OFS="\t" '{print $1,$2,$2+1}' $site_bed | sort -k1,1 -k2,2n | uniq > $temp_bed_1
bedtools intersect -a $temp_bed_1 -b $alu_bed -wo > $temp_bed_2
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_site[$1"\t"$2]++ }ARGIND==2{if($1"\t"$2 in hash_site) hash_snp[$1"\t"$2]++}ARGIND==3{if($1"\t"$2 in hash_site) hash_m6a[$1"\t"$2]++}ARGIND==4{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==5{hash_alu[$1"\t"$2]++}ARGIND==6{site=$1"\t"$2; db="not_in_REDIportal"; snp="non_snp"; m6a="non_m6A_motif"; if(site in hash_db) db="REDIportal"; if(site in hash_snp) snp="snp"; if(site in hash_m6a) m6a="m6A_motif"; label="non_Alu"; if(site in hash_alu) label="Alu"; if($4>=5) print $0,db,snp,m6a,label; }' $site_bed $snp_bed $m6A_motif_bed $REDI_txt $temp_bed_2 $site_bed > $site_tab
if filter_snp=="True" 
then
    mv $site_tab $site_tab.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($8=="non_snp") print $0 }' $site_tab.temp > $site_tab
    rm $site_tab.temp
fi
if filter_m6A=="True" 
then
    mv $site_tab $site_tab.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($9=="non_m6A_motif") print $0 }' $site_tab.temp > $site_tab
    rm $site_tab.temp
fi