site_bed=${1}
alu_bed=${2}
snp_bed=${3}
m6A_motif_bed=${4}
REDI_txt=${5}
temp_bed_1=${6}
temp_bed_2=${7}
site_tab=${8}
cdna2genome_tab=${9}
annotation_gpd=${10}
filter_snp=${11}
filter_m6A=${12}
in_alu_coverage=${13}
out_alu_coverage=${14}
in_alu_ratio=${15}
out_alu_ratio=${16}
coverage_cutoff=${17}
ratio_cutoff=${18}

#add genenome info 10 min
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$1"\t"$2]++}ARGIND==2{if($1"\t"$2 in hash_cdna){ hash_cdna_site[$1"\t"$2]=$3"\t"$4; hash_cdna_strand[$1"\t"$2]=$5;}}ARGIND==3{hash_gene[$2]=$1}ARGIND==4{site=hash_cdna_site[$1"\t"$2]; print $0,site,hash_gene[$1],hash_cdna_strand[$1"\t"$2]; }' $site_bed $cdna2genome_tab $annotation_gpd $site_bed > $site_tab.temp

if test -s $alu_bed; 
then #filter by ALU
   echo "filtering by ALU"
   awk -F '\t' -v OFS="\t" '{if ($8) print $7,$8,$8+1}' $site_tab.temp | sort -k1,1 -k2,2n | uniq > $temp_bed_1
   bedtools intersect -a $temp_bed_1 -b $alu_bed -wo > $temp_bed_2
   
   awk -F '\t' -v OFS="\t" -v in_alu_coverage=$in_alu_coverage -v out_alu_coverage=$out_alu_coverage -v in_alu_ratio=$in_alu_ratio -v out_alu_ratio=$out_alu_ratio 'ARGIND==1{hash_alu[$1"\t"$2]++}ARGIND==2{site=$7"\t"$8; if(site in hash_alu && $4>=in_alu_coverage && $5>=in_alu_ratio) print $0; else if($4>=out_alu_coverage && $5>=out_alu_ratio) print $0; }' $temp_bed_2 $site_tab.temp > $site_tab
else
    cp $alu_bed $temp_bed_2 #an empty file
    cp $alu_bed $temp_bed_1
    awk -F '\t' -v OFS="\t" -v coverage_cutoff=$coverage_cutoff -v ratio_cutoff=$ratio_cutoff 'ARGIND==1{ if($4>=coverage_cutoff && $5>=ratio_cutoff) print $0; }' $site_tab.temp > $site_tab
fi

##annotation on filtered sites

mv $site_tab $site_tab.temp
#head $site_tab.temp

awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$1"\t"$2]++}ARGIND==2{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==3{if($1"\t"$2 in hash_site) hash_snp[$1"\t"$2]++}ARGIND==4{if($1"\t"$2 in hash_site) hash_m6a[$1"\t"$2]++}ARGIND==5{hash_alu[$1"\t"$2]++}ARGIND==6{site=$1"\t"$2; db="not_in_REDIportal"; snp="non_snp"; m6a="non_m6A_motif"; if(site in hash_db) db="REDIportal"; if(site in hash_snp) snp="snp"; label="non_Alu"; if(site in hash_m6a) m6a="m6A_motif"; if(site in hash_alu) label="Alu"; print $0,db,snp,m6a,label; }' $site_tab.temp $REDI_txt $snp_bed $m6A_motif_bed $temp_bed_2 $site_tab.temp  > $site_tab
rm $site_tab.temp

if [ $filter_snp = True ];
then
    echo "filtering out SNPs in [ref_snp_file]"
    mv $site_tab $site_tab.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($12=="non_snp") print $0 }' $site_tab.temp > $site_tab
    rm $site_tab.temp
fi

if [ $filter_m6A = True ];
then
    echo "filtering out m6A motifs generated from [ref_genome_file]"
    mv $site_tab $site_tab.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($13=="non_m6A_motif") print $0 }' $site_tab.temp > $site_tab
    rm $site_tab.temp
fi
