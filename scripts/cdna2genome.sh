cdna_molecule_input=$1
genome_molecule_output=$2
cdna2genome_tab=$3
REDI_txt=$4
annotation_gpd=$5
awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$3"\t"$4]++}ARGIND==2{if($1"\t"$2 in hash_cdna){ hash_cdna_site[$1"\t"$2]=$3"\t"$4; hash_cdna_strand[$1"\t"$2]=$5; hash_site[$3"\t"$4]++}}ARGIND==3{if($1"\t"$2 in hash_site) hash_db[$1"\t"$2]++}ARGIND==4{hash_gene[$2]=$1}ARGIND==5{site=hash_cdna_site[$3"\t"$4]; db="not_in_REDIportal"; if(site in hash_db) db="REDIportal"; print $3,$4,$2,$6,site,hash_gene[$3],hash_cdna_strand[$3"\t"$4],db; }' $cdna_molecule_input $cdna2genome_tab $REDI_txt $annotation_gpd $cdna_molecule_input > $genome_molecule_output