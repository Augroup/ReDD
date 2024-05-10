#!/usr/bin/perl

use strict;
use warnings;

my @sample=("H1-hESC", "DE-H1", "AFG-H1","H9-hESC", "DE-H9", "AFG-H9", "PGC-H1");

foreach my $each (@sample){
	open IN, "$each.ALUAnno.filter.txt" or die $!;
	open OUT, ">$each.meanEditLevel4Isoform";
	print OUT "Isoform\tTotalEditSiteCount\tSampleEditSiteCount\tMeanEditLevel\n";
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
            my $en=0; # total editing sites in this sample
            foreach my $a (@iso){
                $n++;
                $total+=$a;
                if ($a>0) {
                	$en++;
                }
            }
            my $mean=$total/$n;

            print OUT "$key\t$pair{$key}\t$n\t$en\t$mean\n";
    }

    close IN;
    close OUT;

}

## combine edit results of three H1 samples by expression

my %hash;
open IN, "isoform_with_Alu.list" or die $!;
while (<IN>) {
	chomp;
	my @arr=split("\t",$_);
	$hash{$arr[0]}=$arr[1];
}
close IN;


my %h1hesc;
my %h1de;
my %h1pe;

my $tpm_cutoff=1;
open OUT, ">H1hESC_DE_PE.meanEdit4Isoform.compare";
print OUT "Isoform\tGene\thESC_editSiteCount\tMeanEditLevel\tExp\tDE_editSiteCount\tMeanEditLevel\tExp\tPE_editSiteCount\tMeanEditLevel\tExp\tAlutag\n";

open IN, "H1-hESC.meanEditLevel4Isoform" or die $!;
while (<IN>) {
	next if (/^Isoform/);
	chomp;
	my @arr=split("\t",$_);
	my $edit=$arr[3]."\t".$arr[4];
	$h1hesc{$arr[0]}=$edit;
}
close IN;

open IN, "DE-H1.meanEditLevel4Isoform" or die $!;
while (<IN>) {
	next if (/^Isoform/);
	chomp;
	my @arr=split("\t",$_);
	my $edit=$arr[3]."\t".$arr[4];
	$h1de{$arr[0]}=$edit;
}
close IN;

open IN, "AFG-H1.meanEditLevel4Isoform" or die $!;
while (<IN>) {
	next if (/^Isoform/);
	chomp;
	my @arr=split("\t",$_);
	my $edit=$arr[3]."\t".$arr[4];
	$h1pe{$arr[0]}=$edit;
}
close IN;

open IN, "H1hESC_DE_PE.ave.salmon.1TPM" or die $!;
while (<IN>) {
	next if (/^isoform/);
	chomp;
	my @arr=split("\t",$_);
	if ($arr[2]>=$tpm_cutoff and $arr[3]>=$tpm_cutoff and $arr[4]>=$tpm_cutoff) { ## compare if this isoform is expressed in all samples
		if (exists $h1hesc{$arr[0]} and exists $h1de{$arr[0]} and exists $h1pe{$arr[0]}) {
			if (exists $hash{$arr[0]}) {
				print OUT "$arr[0]\t$arr[1]\t$h1hesc{$arr[0]}\t$arr[2]\t$h1de{$arr[0]}\t$arr[3]\t$h1pe{$arr[0]}\t$arr[4]\t$hash{$arr[0]}\n";
			}
			else{
				print OUT "$arr[0]\t$arr[1]\t$h1hesc{$arr[0]}\t$arr[2]\t$h1de{$arr[0]}\t$arr[3]\t$h1pe{$arr[0]}\t$arr[4]\tNonALU\n";
			}
			
		}
	}
}

close IN;
close OUT;


## combine edit results of three H9 samples by expression
my %h9hesc;
my %h9de;
my %h9pe;

open OUT, ">H9hESC_DE_PE.meanEdit4Isoform.compare";
print OUT "Isoform\tGene\thESC_editSiteCount\tMeanEditLevel\tExp\tDE_editSiteCount\tMeanEditLevel\tExp\tPE_editSiteCount\tMeanEditLevel\tExp\tAlutag\n";

open IN, "H9-hESC.meanEditLevel4Isoform" or die $!;
while (<IN>) {
	next if (/^Isoform/);
	chomp;
	my @arr=split("\t",$_);
	my $edit=$arr[3]."\t".$arr[4];
	$h9hesc{$arr[0]}=$edit;
}
close IN;

open IN, "DE-H9.meanEditLevel4Isoform" or die $!;
while (<IN>) {
	next if (/^Isoform/);
	chomp;
	my @arr=split("\t",$_);
	my $edit=$arr[3]."\t".$arr[4];
	$h9de{$arr[0]}=$edit;
}
close IN;

open IN, "AFG-H9.meanEditLevel4Isoform" or die $!;
while (<IN>) {
	next if (/^Isoform/);
	chomp;
	my @arr=split("\t",$_);
	my $edit=$arr[3]."\t".$arr[4];
	$h9pe{$arr[0]}=$edit;
}
close IN;

open IN, "H9hESC_DE_PE.ave.salmon.1TPM" or die $!;
while (<IN>) {
	next if (/^isoform/);
	chomp;
	my @arr=split("\t",$_);
	if ($arr[2]>=$tpm_cutoff and $arr[3]>=$tpm_cutoff and $arr[4]>=$tpm_cutoff) { ## compare if this isoform is expressed in all samples
		if (exists $h9hesc{$arr[0]} and exists $h9de{$arr[0]} and exists $h9pe{$arr[0]}) {
			if (exists $hash{$arr[0]}) {
				print OUT "$arr[0]\t$arr[1]\t$h9hesc{$arr[0]}\t$arr[2]\t$h9de{$arr[0]}\t$arr[3]\t$h9pe{$arr[0]}\t$arr[4]\t$hash{$arr[0]}\n";
			}
			else{
				print OUT "$arr[0]\t$arr[1]\t$h9hesc{$arr[0]}\t$arr[2]\t$h9de{$arr[0]}\t$arr[3]\t$h9pe{$arr[0]}\t$arr[4]\tNonALU\n";
			}
		}
	}
}

close IN;
close OUT;

