mole_txt=$1
site_bed=$2
awk -F '\t' -v OFS="\t" '{hash_total[$3"\t"$4]++; if($6>=0.5) hash_edit[$3"\t"$4]++; hash_score[$3"\t"$4]+=$6; }END{for(site in hash_total){edit=0;if(site in hash_edit) edit=hash_edit[site]; ratio=edit/hash_total[site]; score=hash_score[site]/hash_total[site]; printf site"\t"edit"\t"hash_total[site]"\t%.3f\t%.3f\n",ratio,score; } }' $mole_txt | sort -k1,1 -k2,2n > $site_bed
