ref_genome=$1
motif_bed=$2
bioawk -c fastx '{ print $name, $seq }' $ref_genome > $ref_genome.tab
perl scripts/analysis_extract_m6A_motif.pl $ref_genome.tab $motif_bed
