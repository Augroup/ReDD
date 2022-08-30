site_bed=$1
alu_bed=$2
snp_bed=$3
m6A_motif_bed=$4
REDI_txt=$5
temp_bed_1=$6
temp_bed_2=$7
site_tab=$8
cdna2genome_tab=$9
annotation_gpd=${10}
filter_snp=${11}
filter_m6A=${12}

awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$1"\t"$2]++}ARGIND==2{if($1"\t"$2 in hash_cdna){ hash_cdna_site[$1"\t"$2]=$3"\t"$4; hash_cdna_strand[$1"\t"$2]=$5; hash_site[$3"\t"$4]++}}ARGIND==3{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==4{hash_gene[$2]=$1}ARGIND==5{if($1"\t"$2 in hash_site) hash_snp[$1"\t"$2]++}ARGIND==6{if($1"\t"$2 in hash_site) hash_m6a[$1"\t"$2]++}ARGIND==7{site=hash_cdna_site[$1"\t"$2]; db="not_in_REDIportal"; snp="non_snp"; m6a="non_m6A_motif"; if(site in hash_db) db="REDIportal"; if(site in hash_snp) snp="snp"; if(site in hash_m6a) m6a="m6A_motif"; print $0,site,hash_gene[$1],hash_cdna_strand[$1"\t"$2],db,snp,m6a; }' $site_bed $cdna2genome_tab $REDI_txt $annotation_gpd $snp_bed $m6A_motif_bed $site_bed > $site_tab.temp
awk -F '\t' -v OFS="\t" '{if ($8) print $7,$8,$8+1}' $site_tab.temp | sort -k1,1 -k2,2n | uniq > $temp_bed_1
bedtools intersect -a $temp_bed_1 -b $alu_bed -wo > $temp_bed_2
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_alu[$1"\t"$2]++}ARGIND==2{site=$7"\t"$8; label="non_Alu"; gene=$9; if(site in hash_alu) label="Alu"; if($4>=5) print $0,label; }' $temp_bed_2 $site_tab.temp > $site_tab

if $filter_snp=="True" 
then
    echo $filter_snp
    mv $site_tab $site_tab.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($8=="non_snp") print $0 }' $site_tab.temp > $site_tab
    rm $site_tab.temp
fi
if $filter_m6A=="True" 
then
    echo $filter_m6A
    mv $site_tab $site_tab.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($9=="non_m6A_motif") print $0 }' $site_tab.temp > $site_tab
    rm $site_tab.temp
fi