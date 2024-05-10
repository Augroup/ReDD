#!/usr/bin/perl
use strict;
use warnings;


## this script is used to compared the editing level and sites between gene isoform and TE-derived isoforms
=pod

`cut -f 1-10 Stem_cell_Isoform_output.TE.annotation.list2 > 001`;
`cut -f 1-10 Stem_cell_Isoform_output.TE-Gene.annotation.list2 | sed '1d' > 002`;
`cat 001 002 > Stem_cell.TE-derived.list`;
`rm 001 002`;
=cut

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

	open SUM, ">$name.TE_vs_gene.cmp.4R";
	print SUM "ID\teditSite\teditLevel\tclass\n";

	open IN, "Stem_cell_Isoform_output.Gene.annotation.list2" or die $!;
	open OUT, ">$name.gene.editsite_editlevel.txt";
	print OUT "isoform\tgene\teditSite\teditLevel\n";
	while (<IN>) {
		next if (/^isoform/);
		my @arr=split("\t",$_);
		if (exists $count{$arr[0]}) {
			print OUT "$arr[0]\t$arr[1]\t$count{$arr[0]}\t$level{$arr[0]}\n";
			print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\tgene\n";
		}
	}
	close IN;
	close OUT;

	open IN, "Stem_cell.TE-derived.list" or die $!;
	open OUT, ">$name.te.editsite_editlevel.txt";
	print OUT "isoform\tgene\teditSite\teditLevel\tTEproportion\tTEsubfamily\tTEfamily\tTEtype\n";
	while (<IN>) {
		next if (/^isoform/);
		chomp;
		my @arr=split("\t",$_);
		if (exists $count{$arr[0]}) {
			print OUT "$arr[0]\t$arr[1]\t$count{$arr[0]}\t$level{$arr[0]}\t$arr[7]\t$arr[8]\t$arr[9]\n";
			print SUM "$arr[0]\t$count{$arr[0]}\t$level{$arr[0]}\tTE\n";
		}
	}
	close IN;
	close OUT;
	close SUM;
}