#!/usr/bin/perl
use strict;
use warnings;

## This script is used to filter out the low-confidence editing sites.

## filtering parameters

my $cov_cutoff=10;  # coverage 
my $ratio_cutoff=0.1; # editing level
my $read_cutoff=2; # reads to sopport the editing status

my $cov_alu=10;  # coverage for alu region
my $ratio_alu=0.1; # editing level for alu region
my $read_alu=2; # reads to sopport the editing status for alu region


my @files=glob("*.ALUAnno");
foreach my $each (@files){
	open IN, "$each" or die $!;
	open OUT, ">$each.filter.txt";
	while (<IN>) {
		chomp;
		my $line=$_;
		my @arr=split("\t",$line);

		if ($arr[-1] eq "NonALU") {
			if ($arr[2]>=$read_cutoff and $arr[3]>=$cov_cutoff and $arr[4]>=$ratio_cutoff) {
				print OUT "$line\n";
			}
		}
		elsif ($arr[-1] eq "ALU") {
			if ($arr[2]>=$read_alu and $arr[3]>=$cov_alu and $arr[4]>=$ratio_alu) {
				print OUT "$line\n";
			}	
		}
	}
	close IN;
	close OUT;
}

my @file2=glob("*.txt");
open OUT, ">ALU_nonALU.sum";
print OUT "Name\tALU_count\tnonALU_count\n";
foreach my $each (@file2){
	open IN, "$each" or die $!;
	my $db=0;
	my $ndb=0;
	while (<IN>) {
		chomp;
		my @arr=split("\t",$_);
		if ($arr[-1] eq "ALU") {
			$db++;
		}
		else{
			$ndb++;
		}
	}
	print OUT "$each\t$db\t$ndb\n";
	close IN;
}
close OUT;

