
awk -F '\t' -v OFS="\t" 'ARGIND==1{ if($8=="non_snp") print $0 }' $genome_output.tab > $genome_output.flt_snp.tab