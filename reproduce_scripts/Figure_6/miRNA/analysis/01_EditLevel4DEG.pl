#!/usr/bin/perl
use strict;
use warnings;

# extract mean editing level for DEGs

my $logfc=1;
my $padj=0.05;


my @files=glob("Sig_DEGs.*.isoform.tab");
foreach my $each (@files){
	$each=~/Sig_DEGs.(.*)_vs_(.*).isoform.tab/;
	my $sample1=$1;
	my $sample2=$2;
	my %hash1;
	my %hash2;
	my %gene;
	open IN, "$sample1.meanEditLevel4Isoform" or die $!;
	while (<IN>) {
		next if (/^Isoform/);
		chomp;
		my @arr=split("\t",$_);
		$hash1{$arr[0]}=$arr[4];
		$gene{$arr[0]}=$arr[1];
	}
	close IN;

	open IN, "$sample2.meanEditLevel4Isoform" or die $!;
	while (<IN>) {
		next if (/^Isoform/);
		chomp;
		my @arr=split("\t",$_);
		$hash2{$arr[0]}=$arr[4];
		$gene{$arr[0]}=$arr[1];
	}
	close IN;

	open IN, "$each" or die $!;
	open OUT, ">$sample1.$sample2.DEG.MeanEditLevel.txt";
	print OUT "Isoform\tgene\tlogFC\tpvalue\tpadj\tMean1\tMean2\tImean1\tImean2\n";
	while (<IN>) {
		next if (/^GeneName/);
		my @arr=split("\t",$_);
		if (abs($arr[2])>=$logfc and $arr[6]<=$padj) {
			if (exists $hash1{$arr[0]} and exists $hash2{$arr[0]}) {
				my $mean1=($arr[7]+$arr[8]+$arr[9])/3;
				my $mean2=($arr[10]+$arr[11]+$arr[12])/3;
				print OUT "$arr[0]\t$gene{$arr[0]}\t$arr[2]\t$arr[5]\t$arr[6]\t$mean1\t$mean2\t$hash1{$arr[0]}\t$hash2{$arr[0]}\n";
			}
			elsif (exists $hash1{$arr[0]} and not exists $hash2{$arr[0]}) {
				my $mean1=($arr[7]+$arr[8]+$arr[9])/3;
				my $mean2=($arr[10]+$arr[11]+$arr[12])/3;
				print OUT "$arr[0]\t$gene{$arr[0]}\t$arr[2]\t$arr[5]\t$arr[6]\t$mean1\t$mean2\t$hash1{$arr[0]}\t0\n";
			}
			elsif (not exists $hash1{$arr[0]} and exists $hash2{$arr[0]}) {
				my $mean1=($arr[7]+$arr[8]+$arr[9])/3;
				my $mean2=($arr[10]+$arr[11]+$arr[12])/3;
				print OUT "$arr[0]\t$gene{$arr[0]}\t$arr[2]\t$arr[5]\t$arr[6]\t$mean1\t$mean2\t0\t$hash2{$arr[0]}\n";
			}
		}
	}
	close IN;
	close OUT;

}
