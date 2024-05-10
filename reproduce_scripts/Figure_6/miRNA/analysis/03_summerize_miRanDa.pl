#!/usr/bin/perl
use strict;
use warnings;

use lib "/home/libowyj/perl5/lib/perl5";
use Array::Utils qw(:all);


my @files=glob("*.fasta.miranda");
foreach my $each (@files){
	open IN, "$each" or die $!;
	open OUT, ">$each.tmp";
	while (<IN>) {
		if (/^\>\>/) {
			print OUT "$_";
		}
	}
	close IN;
	close OUT;
}

my $score=145;
my $energy=-20;

my @sample=("H1-hESC", "DE-H1", "AFG-H1","H9-hESC", "DE-H9", "AFG-H9");

foreach my $each (@sample){
	my %before_all;
	my %before_above;
	my %before_below;
	open IN, "$each.before.fasta.miranda.tmp" or die $!;
	while (<IN>) {
		s/\>\>//;
		my @arr=split("\t",$_);
		push @{$before_all{$arr[1]}},$arr[0];
		if ($arr[2]>=$score and $arr[3]<=$energy) {
			push @{$before_above{$arr[1]}},$arr[0];
		}
		else{
			push @{$before_below{$arr[1]}},$arr[0];
		}
	}
	close IN;

	my %after_all;
	my %after_above;
	my %after_below;
	open IN, "$each.after.fasta.miranda.tmp" or die $!;
	while (<IN>) {
		s/\>\>//;
		my @arr=split("\t",$_);
		push @{$after_all{$arr[1]}},$arr[0];
		if ($arr[2]>=$score and $arr[3]<=$energy) {
			push @{$after_above{$arr[1]}},$arr[0];
		}
		else{
			push @{$after_below{$arr[1]}},$arr[0];
		}
	}
	close IN;

	open OUT, ">$each.miRNA_changes_before_and_after_editing.gain_loss.sum";
	foreach my $key (keys %after_above){
		my @list_after=@{$after_above{$key}};
		if (exists $before_all{$key}) {
			my @list_before=@{$before_all{$key}};
			my @minus1 = array_minus( @list_after, @list_before);
			my $minus1=@minus1;
			if ($minus1>0) {
				my $tmp=join("\t",@minus1);
				print OUT "$key\t$tmp\tmiRNA_gain\n";
			}
		}
	}

	foreach my $key (keys %before_above){
		my @list_before=@{$before_above{$key}};
		if (exists $after_all{$key}) {
			my @list_after=@{$after_all{$key}};
			my @minus2 = array_minus( @list_before, @list_after);
			my $minus2=@minus2;
			if ($minus2>0) {
				my $tmp=join("\t",@minus2);
				print OUT "$key\t$tmp\tmiRNA_loss\n";
			}
		}
	}
	close OUT;


	open OUT, ">$each.miRNA_changes_before_and_after_editing.strengthen_weaken.sum";
	foreach my $key (keys %after_above){
		my @list_after=@{$after_above{$key}};
		if (exists $before_below{$key}) {
			my @list_before=@{$before_below{$key}};
			my @minus1 = intersect( @list_after, @list_before);
			my $minus1=@minus1;
			if ($minus1>0) {
				my $tmp=join("\t",@minus1);
				print OUT "$key\t$tmp\tmiRNA_up\n";
			}
		}
	}

	foreach my $key (keys %before_above){
		my @list_before=@{$before_above{$key}};
		if (exists $after_below{$key}) {
			my @list_after=@{$after_below{$key}};
			my @minus1 = intersect(@list_before, @list_after);
			my $minus1=@minus1;
			if ($minus1>0) {
				my $tmp=join("\t",@minus1);
				print OUT "$key\t$tmp\tmiRNA_down\n";
			}
		}
	}

	close OUT;
}


## do summary
open OUT, ">StemCell.miRNA_binding.gain_loss.sum" or die $!;
print OUT "cell\tisoformCount\tisoformCountGain\ttotalGain\tmeanGain\tisoformCountLoss\ttotalLoss\tmeanLoss\n";
open OUT2, ">StemCell.miRNA_binding.gain_loss.4R" or die $!;
print OUT2 "cell\tisoform\tcount\ttype\n";
my @out1=glob("*.miRNA_changes_before_and_after_editing.gain_loss.sum");
foreach my $each (@out1){
	$each=~/(.*).miRNA_changes_before_and_after_editing.gain_loss.sum/;
	my $cell=$1;
	my $isoform_count=0;
	my $isoform_gain=0;
	my $isoform_loss=0;
	my $gain=0;
	my $loss=0;
	open IN, "$each" or die $!;
	while (<IN>) {
		chomp;
		$isoform_count++;
		my @arr=split("\t",$_);
		my $id=shift @arr;
		my $tag=pop @arr;
		my $arr=@arr;
		if ($tag eq "miRNA_gain") {
			print OUT2 "$cell\t$id\t$arr\tmiRNA_gain\n";
			$gain+=$arr;
			$isoform_gain++;
		}
		elsif ($tag eq "miRNA_loss") {
			print OUT2 "$cell\t$id\t$arr\tmiRNA_loss\n";
			$loss+=$arr;
			$isoform_loss++;
		}
	}
	my $mean_gain=$gain/$isoform_gain;
	my $mean_loss=$loss/$isoform_loss;
	print OUT "$cell\t$isoform_count\t$isoform_gain\t$gain\t$mean_gain\t$isoform_loss\t$loss\t$mean_loss\n";
} 

close OUT;
close OUT2;


open OUT, ">StemCell.miRNA_binding.strengthen_weaken.sum" or die $!;
print OUT "cell\tisoformCount\tisoformCountGain\ttotalGain\tmeanGain\tisoformCountLoss\ttotalLoss\tmeanLoss\n";
open OUT2, ">StemCell.miRNA_binding.strengthen_weaken.4R" or die $!;
print OUT2 "cell\tisoform\tcount\ttype\n";
my @out2=glob("*.miRNA_changes_before_and_after_editing.strengthen_weaken.sum");
foreach my $each (@out2){
	$each=~/(.*).miRNA_changes_before_and_after_editing.strengthen_weaken.sum/;
	my $cell=$1;
	my $isoform_count=0;
	my $isoform_gain=0;
	my $isoform_loss=0;
	my $gain=0;
	my $loss=0;
	open IN, "$each" or die $!;
	while (<IN>) {
		chomp;
		$isoform_count++;
		my @arr=split("\t",$_);
		my $id=shift @arr;
		my $tag=pop @arr;
		my $arr=@arr;
		if ($tag eq "miRNA_up") {
			print OUT2 "$cell\t$id\t$arr\tmiRNA_up\n";
			$gain+=$arr;
			$isoform_gain++;
		}
		elsif ($tag eq "miRNA_down") {
			print OUT2 "$cell\t$id\t$arr\tmiRNA_down\n";
			$loss+=$arr;
			$isoform_loss++;
		}
	}
	my $mean_gain=$gain/$isoform_gain;
	my $mean_loss=$loss/$isoform_loss;
	print OUT "$cell\t$isoform_count\t$isoform_gain\t$gain\t$mean_gain\t$isoform_loss\t$loss\t$mean_loss\n";
} 

close OUT;
close OUT2;




