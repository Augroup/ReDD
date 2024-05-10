#!/usr/bin/perl
use strict;
use warnings;

open SUM, ">TE_family.trans_count.sum";

my @files=glob("*.meanEditLevel4Isoform");
foreach my $each (@files){
	my %count;
	my %level;
	$each=~/(.*).meanEditLevel4Isoform/;
	my $name=$1;
	open IN, "$each" or die $!;
	while (<IN>) {
		next if (/^Isoform/);
		chomp;
		my @arr=split("\t",$_);
		$count{$arr[0]}=$arr[2];
		$level{$arr[0]}=$arr[3];
	}
	close IN;

	open OUT, ">$name.gene.TEfamily.4test";
	#print OUT "gene\n";
	open IN, "Stem_cell_Isoform_output.Gene.annotation.list2" or die $!;
	while (<IN>) {
		next if (/^isoform/);
		my @arr=split("\t",$_);
		if (exists $count{$arr[0]}) {
			print OUT "$level{$arr[0]}\n";
		}
	}
	close IN;

	open IN, "Stem_cell.TE-derived.list" or die $!;
	my %fam;
	while (<IN>) {
		next if (/^isoform/);
		chomp;
		my @arr=split("\t",$_);
		my $te=$arr[8]."_".$arr[9];
		push @{$fam{$te}}, $arr[0];
		
	}

	foreach my $key (keys %fam){
		my @each=@{$fam{$key}};
		my $each=@each;
		next if ($each<50); 
		open OUT2, ">$name.$key.TEfamily.4test";
		#print OUT2 "$key\n";
		my $n=0; ## the number of transcript for each family
		foreach my $a (@each){
			if (exists $level{$a}) {
				$n++;
				print OUT2 "$level{$a}\n";
			}
		}
		close OUT2;
		print SUM "$name\t$key\t$n\n";
	}
	close IN;
	close OUT;
}

close SUM;

`mkdir -p family_analysis`;
`mv *.TEfamily.4test family_analysis`;