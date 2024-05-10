#!/usr/bin/perl

use strict;
use warnings;

my @sample=("H1-hESC", "DE-H1", "AFG-H1","H9-hESC", "DE-H9", "AFG-H9", "PGC-H1");

foreach my $each (@sample){
	open IN, "$each.meanEditLevel4Isoform" or die $!;
	open OUT, ">$each.editDiversity";
	print OUT "gene\tcv\tfc\n";
	my %gene;
	while (<IN>) {
		next if (/^Isoform/);
		chomp;
		my @arr=split("\t",$_);
		push @{$gene{$arr[1]}},$arr[4];
	}

	foreach my $key (keys %gene){
		my @value=@{$gene{$key}};
		my $value=@value;
		my @sort=sort @value;
		next if $sort[0]==0;
		my $fc=$sort[-1]/$sort[0];
		next if $value==1; ## next if single isoform genes
		my $ave=average(\@value);
		my $std=stdev(\@value);
		if ($ave==0) {
			next;
			#print OUT "0\t$pop[1]\n";
			#print OUT3 "$gene\t0\n";
		}
		else{
			my $cv=$std/$ave;
			print OUT "$key\t$cv\t$fc\n";
		}	
	}
	close IN;
	close OUT;
}

sub average {
@_ == 1 or die ('Sub usage: $average = average(\@array);');
my ($array_ref) = @_;
my $sum;
my $count = scalar @$array_ref;
foreach (@$array_ref) { $sum += $_; }
return $sum / $count;
}

sub stdev{
        @_ == 1 or die ('Sub usage: $average = average(\@array);');
		my ($array_ref) = @_;
        my $average = &average($array_ref);
        my $sqtotal = 0;
        foreach(@$array_ref) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$array_ref-1)) ** 0.5;
        return $std;
}