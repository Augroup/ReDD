site_bed=${1}
alu_bed=${2}
snp_bed=${3}
m6A_motif_bed=${4}
REDI_txt=${5}
temp_bed_1=${6}
temp_bed_2=${7}
site_tab=${8}
filter_snp=${9}
filter_m6A=${10}
in_alu_coverage=${11}
out_alu_coverage=${12}
in_alu_ratio=${13}
out_alu_ratio=${14}
coverage_cutoff=${15}
ratio_cutoff=${16}


if test -s $alu_bed; 
then #filter by ALU
   echo "filtering by ALU in [ref_alu]"
   awk -F '\t' -v OFS="\t" '{print $1,$2,$2+1}' $site_bed | sort -k1,1 -k2,2n | uniq > $temp_bed_1
   bedtools intersect -a $temp_bed_1 -b $alu_bed -wo > $temp_bed_2
   
   awk -F '\t' -v OFS="\t" -v in_alu_coverage=$in_alu_coverage -v out_alu_coverage=$out_alu_coverage -v in_alu_ratio=$in_alu_ratio -v out_alu_ratio=$out_alu_ratio 'ARGIND==1{hash_alu[$1"\t"$2]++}ARGIND==2{site=$1"\t"$2; if(site in hash_alu && $4>=in_alu_coverage && $5>=in_alu_ratio) print $0; else if($4>=out_alu_coverage && $5>=out_alu_ratio) print $0; }' $temp_bed_2 $site_bed  > $site_tab

else #filter by coverage_cutoff and ratio_cutoff
    cp $alu_bed $temp_bed_2 #an empty file
    cp $alu_bed $temp_bed_1
    awk -F '\t' -v OFS="\t" -v coverage_cutoff=$coverage_cutoff -v ratio_cutoff=$ratio_cutoff 'ARGIND==1{ if($4>=coverage_cutoff && $5>=ratio_cutoff) print $0; }' $site_bed > $site_tab
fi

#annotation on filtered sites
mv $site_tab $site_tab.temp
#head $site_tab.temp
#for debug
#echo $site_tab.temp 
#echo $snp_bed 
#echo $m6A_motif_bed 
#echo $REDI_txt 
#echo $temp_bed_2

awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_site[$1"\t"$2]++ }ARGIND==2{if($1"\t"$2 in hash_site) hash_snp[$1"\t"$2]++}ARGIND==3{if($1"\t"$2 in hash_site) hash_m6a[$1"\t"$2]++}ARGIND==4{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==5{hash_alu[$1"\t"$2]++}ARGIND==6{site=$1"\t"$2; db="not_in_REDIportal"; snp="non_snp"; m6a="non_m6A_motif"; if(site in hash_db) db="REDIportal"; if(site in hash_snp) snp="snp"; if(site in hash_m6a) m6a="m6A_motif"; label="non_Alu"; if(site in hash_alu) label="Alu"; print $0,db,snp,m6a,label;}' $site_tab.temp $snp_bed $m6A_motif_bed $REDI_txt $temp_bed_2 $site_tab.temp > $site_tab
#head $site_tab
rm $site_tab.temp
if [ $filter_snp = True ];
then
    echo "filtering out SNPs in [ref_snp_file]"
    mv $site_tab $site_tab.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($8=="non_snp") print $0 }' $site_tab.temp > $site_tab
    rm $site_tab.temp
fi
if [ $filter_m6A = True ];
then
    echo "filtering out m6A motifs generated from [ref_genome_file]"
    mv $site_tab $site_tab.temp
    awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($9=="non_m6A_motif") print $0 }' $site_tab.temp > $site_tab
    rm $site_tab.temp
fi
