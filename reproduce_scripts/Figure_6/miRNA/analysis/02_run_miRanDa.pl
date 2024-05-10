#!/usr/bin/perl
use strict;
use warnings;


## prepare miRNA and mRNA sequence
my %miRNA;
my %isoform;
$/=">";
open IN, "hsa.miRNA.fa" or die $!;
while (<IN>) {
	next if (length $_<2);
	s/>//;
	my @arr=split("\n",$_);
	$miRNA{$arr[0]}=$arr[1];
}
close IN;

open IN, "Stem_cell_talon.flt.bam_flt.gtf.fa" or die $!;
while (<IN>) {
	next if (length $_<2);
	s/>//;
	my $line=$_;
	my @arr=split("\n",$line);
	$isoform{$arr[0]}=$line;
}
close IN;

$/="\n";

open IN, "hESC_DE.expressed.miRNA.list" or die $!;
open OUT, ">stemcell.miRNA.fa";
while (<IN>) {
	chomp;
	my @arr=split("\t",$_);
	next if ($arr[1]<=100 and $arr[2]<=100); # keep highly expressed miRNA
	if (exists $miRNA{$arr[0]}) {
		print OUT ">$arr[0]\n$miRNA{$arr[0]}\n";
	}
}
close IN;
close OUT;

my @files=glob("*.MeanEditLevel.txt");
foreach my $each (@files){
	$each=~/(.*).DEG.MeanEditLevel.txt/;
	my $sample=$1;
	open IN, "$each" or die $!;
	open OUT, ">$sample.isoform.fa";
	while (<IN>) {
		next if (/^Isoform/);
		my @arr=split("\t",$_);
		if (exists $isoform{$arr[0]}) {
			print OUT ">$isoform{$arr[0]}";
		}
	}
	close IN;
	close OUT;

}

## convert A to G for editing sites on mRNA sequences

my @fa=glob("*.isoform.fa");
foreach my $each (@fa){
	$each=~/(.*)\.(.*)\.isoform\.fa/;
	my $sample1=$1;
	my $sample2=$2;

	my %pos1;
	open IN, "$sample1.ALUAnno.filter.txt" or die $!;
	while (<IN>) {
		my @arr=split("\t",$_);
		push @{$pos1{$arr[0]}},$arr[1];
	}
	close IN;

	my %pos2;
	open IN, "$sample2.ALUAnno.filter.txt" or die $!;
	while (<IN>) {
		my @arr=split("\t",$_);
		push @{$pos2{$arr[0]}},$arr[1];
	}
	close IN;

	$/=">";
	open IN, "$each" or die $!;
	open OUT, ">$sample1.before.fasta";
	open OUT2, ">$sample1.after.fasta";
	open OUT3, ">$sample2.before.fasta";
	open OUT4, ">$sample2.after.fasta";
	while (<IN>) {
		next if (length $_<2);
		s/>//;
		my @arr=split("\n",$_);
		my $id=shift @arr;
		my $seq=join("",@arr);
		my $new=$seq;
		if (exists $pos1{$id}) {
			my @list1=@{$pos1{$id}};
			foreach my $a (@list1){
				my $site=$a-1;
				my $new=substr($new, $site, 1, "G"); ## replace the editing sites from A into G 
			}
			print OUT ">$id\n$seq\n";
			print OUT2 ">$id\n$new\n";
		}
		$new=$seq;
		if (exists $pos2{$id}) {
			my @list2=@{$pos2{$id}};
			foreach my $a (@list2){
				my $site=$a-1;
				my $new=substr($new, $site, 1, "G"); ## replace the editing sites from A into G 
			}
			print OUT3 ">$id\n$seq\n";
			print OUT4 ">$id\n$new\n";
		}
	}
	close IN;
	close OUT;
	close OUT2;
	close OUT3;
	close OUT4;
	$/="\n";
}

my @file2=glob("*.fasta");
foreach my $each (@file2){
	`~/local/miRanda-1.9/bin/miranda stemcell.miRNA.fa $each -out $each.miranda`;
}







