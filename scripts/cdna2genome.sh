cdna_molecule_input=$1
genome_molecule_output=$2
genome_site_output=$3
cdna2genome_tab=$4
cdna_site_tab=$5

awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_cdna[$3"\t"$4]++}ARGIND==2{if($1"\t"$2 in hash_cdna){ hash_cdna_site[$1"\t"$2]=$3"\t"$4; hash_cdna_strand[$1"\t"$2]=$5; hash_site[$3"\t"$4]++}}ARGIND==3{site=hash_cdna_site[$3"\t"$4]; print "I",$2,$3,$4,hash_cdna_strand[$3"\t"$4],$6,site; }' $cdna_molecule_input $cdna2genome_tab $cdna_molecule_input > $genome_molecule_output

awk -F '\t' -v OFS="\t" 'ARGIND==1{hash_site[$8"\t"$9]=$11"\t"$12"\t"$13"\t"$14}ARGIND==2{hash_total[$7"\t"$8]++; if($6>=0.5) hash_edit[$7"\t"$8]++; hash_score[$7"\t"$8]+=$6; }END{for(site in hash_total){if(site in hash_site) {edit=0;if(site in hash_edit) edit=hash_edit[site]; ratio=edit/hash_total[site]; score=hash_score[site]/hash_total[site]; print site,edit,hash_total[site],ratio,score,hash_site[site]}}}' $cdna_site_tab $genome_molecule_output | sort -k1,1 -k2,2n > $genome_site_output

