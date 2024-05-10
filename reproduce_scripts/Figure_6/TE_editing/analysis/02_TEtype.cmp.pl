#!/usr/bin/perl
use strict;
use warnings;

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

	open SUM, ">$name.TE_vs_gene.TEtype.4R";
	print SUM "ID\teditSite\teditLevel\tclass\n";

	open IN, "Stem_cell_Isoform_output.Gene.annotation.list2" or die $!;
	while (<IN>) {
		next if (/^isoform/);
		my @arr=split("\t",$_);
		if (exists $count{$arr[0]}) {
			print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\tgene\n";
		}
	}
	close IN;

	open IN, "Stem_cell.TE-derived.list" or die $!;
	while (<IN>) {
		next if (/^isoform/);
		chomp;
		my @arr=split("\t",$_);
		if (exists $count{$arr[0]}) {
			if ($arr[-1] eq "DNA") {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\tDNA\n";
			}
			elsif ($arr[-1] eq "LINE") {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\tLINE\n";
			}
			elsif ($arr[-1] eq "SINE") {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\tSINE\n";
			}
			elsif ($arr[-1] eq "LTR") {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\tLTR\n";
			}
			elsif ($arr[-1] eq "Retroposon") {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\tRetroposon\n";
			}
		}
	}
	close IN;
	close SUM;
}