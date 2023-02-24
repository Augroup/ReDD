#!/usr/bin/perl -w
##usage: perl gpd cdna2genome.tab
##tranfer gpd file to cdna2genome file
$gpd=$ARGV[0];
$cdna2genome=$ARGV[1];

##cdna2genome
open(IN, $gpd) or die("Could not open file IN $gpd.\n");
open(OUT, "> $cdna2genome") or die("Could not open file $cdna2genome.\n");
while($line=<IN>){
	chomp $line;
	@data=split /\t/,$line;
	$exon_num=$data[8];
	$exon_l = $data[9];
	$exon_r = $data[10];
	@exon_l=split /,/,$exon_l;
	@exon_r=split /,/,$exon_r;
	$iso_len=0;
	for($i=1;$i<=$exon_num;$i++){
		$iso_len=$iso_len+($exon_r[($i-1)]-$exon_l[($i-1)]);
	}
	$cdna_len=0;
	for($i=1;$i<=$exon_num;$i++){
		for($j=$exon_l[($i-1)]+1;$j<$exon_r[($i-1)]+1;$j++){
			$cdna_len++;
			if($data[3] eq "+"){
				$seq=$cdna_len;
			}else{
				$seq=$iso_len+1-$cdna_len;
			}
			print OUT $data[1]."\t".$seq."\t".$data[2]."\t".$j."\t".$data[3]."\n";
		}
	}
}
close IN;
close OUT;
print "Done\n";