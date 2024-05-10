#!/usr/bin/perl
use strict;
use warnings;

## This script is used to analyze the RNA editing profile with TE-Gene transcripts

my @allfiles=glob("*.ALUAnno.filter.txt");
foreach my $file (@allfiles){
	$file=~/(.*).ALUAnno.filter.txt/;
	my $name=$1;
	my %hash;
	open IN, "$file" or die $!;
	while (<IN>) {
		chomp;
		my @arr=split("\t",$_);
		my $id=$arr[0]."_".$arr[1];
		$hash{$id}=$arr[4];
	}
	close IN;

	open IN, "TE-Gene.TE_nonTE.bed" or die $!;
	open OUT, ">$name.TE-Gene.TE_nonTE.editLevel.txt";
	print OUT "Isoform\tStart\tEnd\tSubfamily\tFamily\tType\tEditCount\tEditLevel\n";
	while (<IN>) {
		chomp;
		next if (/LTR\?/);
		my $line=$_;
		my @arr=split("\t",$line);
		my $edit_a=0;
		my $total_ratio=0;
		my $mean=0;
		for (my $n = $arr[1]; $n < $arr[2]; $n++) {
			my $site=$arr[0]."_".$n;
			if (exists $hash{$site}) {
				$edit_a++;
				$total_ratio+=$hash{$site};
			}
		}
		if ($edit_a>0) {
			$mean=$total_ratio/$edit_a;
			print OUT "$line\t$edit_a\t$mean\n";
		}
		else{
			$mean=0;
			print OUT "$line\t$edit_a\t$mean\n";
		}
	}
	close IN;
	close OUT;
}



sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}






