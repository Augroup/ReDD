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

	open SUM, ">$name.TE_vs_gene.proprotion.4R";
	print SUM "ID\teditSite\teditLevel\tclass\n";

	open IN, "Stem_cell_Isoform_output.Gene.annotation.list2" or die $!;
	while (<IN>) {
		next if (/^isoform/);
		my @arr=split("\t",$_);
		if (exists $count{$arr[0]}) {
			print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t0\n";
		}
	}
	close IN;

	open IN, "Stem_cell.TE-derived.list" or die $!;
	while (<IN>) {
		next if (/^isoform/);
		chomp;
		my @arr=split("\t",$_);
		if (exists $count{$arr[0]}) {
			if ($arr[4]<=0.1) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t10\n";
			}
			elsif ($arr[4]>0.1 and $arr[4]<=0.2) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t20\n";
			}
			elsif ($arr[4]>0.2 and $arr[4]<=0.3) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t30\n";
			}
			elsif ($arr[4]>0.3 and $arr[4]<=0.4) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t40\n";
			}
			elsif ($arr[4]>0.4 and $arr[4]<=0.5) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t50\n";
			}
			elsif ($arr[4]>0.5 and $arr[4]<=0.6) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t60\n";
			}
			elsif ($arr[4]>0.6 and $arr[4]<=0.7) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t70\n";
			}
			elsif ($arr[4]>0.7 and $arr[4]<=0.8) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t80\n";
			}
			elsif ($arr[4]>0.8 and $arr[4]<=0.9) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t90\n";
			}
			elsif ($arr[4]>0.9 and $arr[4]<=1) {
				print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\t100\n";
			}
			
		}
	}
	close IN;
	close SUM;
}