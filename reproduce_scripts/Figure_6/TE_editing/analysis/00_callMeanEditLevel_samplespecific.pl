#!/usr/bin/perl

use strict;
use warnings;

my @sample=("H1-hESC", "DE-H1", "AFG-H1","H9-hESC", "DE-H9", "AFG-H9", "PGC-H1");

foreach my $each (@sample){
	open IN, "$each.ALUAnno.filter.txt" or die $!;
	open OUT, ">$each.meanEditLevel4Isoform";
	print OUT "Isoform\tgene\tSampleEditSiteCount\tMeanEditLevel\n";
	my %pair; # isoform-gene pair
	my %isoform;
	while (<IN>) {
		chomp;
		my @arr=split("\t",$_);
		$pair{$arr[0]}=$arr[8];
		push @{$isoform{$arr[0]}}, $arr[4];
	}

	foreach my $key (keys %isoform){
            my @iso=@{$isoform{$key}};
            my $total=0;
            my $n=0; # total editing sites on this isoform
            foreach my $a (@iso){
                $n++;
                $total+=$a;
            }
            my $mean=$total/$n;

            print OUT "$key\t$pair{$key}\t$n\t$mean\n";
    }

    close IN;
    close OUT;
}
