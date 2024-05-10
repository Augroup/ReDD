#!/usr/bin/perl
use strict;
use warnings;

## This script is to compare the editing level among single copy, multiple copy and sense-antisense paired Alu copies.
open OUT, ">Alu_type.list";
print OUT "ID\tType\n";

my %hash;
open IN, "trans_with_sense_antisense.Alu.out" or die $!;  
while (<IN>) {
	my @arr=split("\t",$_);
	$hash{$arr[0]}=1;
}
close IN;

foreach my $key (keys %hash){
	print OUT "$key\tPaired\n";
}

my %hash2;
open IN, "Alu.multipleCopy.list" or die $!;
while (<IN>) {
	my @arr=split("\t",$_);
	$hash2{$arr[0]}=1;
	if (not exists $hash{$arr[0]}) {
		print OUT "$arr[0]\tMulti\n";
	}
}
close IN;

open IN, "Alu.singleCopy.list" or die $!;
while (<IN>) {
	my @arr=split("\t",$_);
	if (not exists $hash{$arr[0]} and not exists $hash2{$arr[0]}) {
		print OUT "$arr[0]\tSingle\n";
	}
}
close IN;
close OUT;


my @files=glob("*.meanEditLevel4Isoform");
foreach my $each (@files){
	open IN, "$each" or die $!;
	$each=~/(.*).meanEditLevel4Isoform/;
	my $sample=$1;
	open OUT, ">$sample.Alu_type.Edit.cmp";
	print OUT "Isoform\tGene\tEditSite\tEditLevel\tType\n";
	my %edit;
	while (<IN>) {
		next if (/^Isoform/);
		chomp;
		my $line=$_;
		my @arr=split("\t",$line);
		$edit{$arr[0]}=$line;
	}
	close IN;

	open IN, "Alu_type.list" or die $!;
	while (<IN>) {
		next if (/^ID/);
		chomp;
		my @arr=split("\t",$_);
		if (exists $edit{$arr[0]}) {
			print OUT "$edit{$arr[0]}\t$arr[1]\n";
		}
	}
	close IN;
	close OUT;
}



