#!/usr/bin/perl -w

# $result_dir="/fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/1-Analysis/threshold_evaluation/hg19tohg38";
no warnings;
$out=$ARGV[1];
$ref=$ARGV[0];
open(IN, $ref) or die("Could not open file IN.\n");
open(OUT, ">".$out) or die("Could not open file OUT.\n");
@motifs=("GGACA","AGACA","AAACA","GAACA","GGACC","AGACC","AAACC","GAACC","GGACT","AGACT","AAACT","GAACT");
while($line=<IN>){
	chomp $line;
	@data=split /\t/,$line;
	$seq=$data[1];
	foreach $motif(@motifs){
		$offset=0;
		while ($pos != -1){
			if ($offset>0){
				print OUT $data[0]."\t".($pos+3)."\t".($pos+4)."\t"."A"."\t".$motif."\n";
			}
			$offset = $pos + 1;
			$pos = index($seq, $motif, $offset);
		}
		$revmotif=reverse $motif;
		$revmotif=~tr/ATGCatgc/TACGtacg/;
		$offset=0;
		$pos1=0;
		while ($pos1 != -1){
			if ($offset>0){
				print OUT $data[0]."\t".($pos1+3)."\t".($pos1+4)."\t"."T"."\t".$revmotif."\n";
			}
			$offset = $pos1 + 1;
			$pos1 = index($seq, $revmotif, $offset);
		}
	}
	print "$data[0] Done\n";
}
close IN;
close OUT;
print "Done\n";