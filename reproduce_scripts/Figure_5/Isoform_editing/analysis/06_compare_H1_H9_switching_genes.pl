#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

my (@h1hesc2de, @h9hesc2de, @h1de2pe, @h9de2pe);

open IN, "H1hESC_DE.editSwitch.txt.4R" or die $!;
while (<IN>) {
	my @arr=split("\t",$_);
	push @h1hesc2de, $arr[0];
}
close IN;

open IN, "H9hESC_DE.editSwitch.txt.4R" or die $!;
while (<IN>) {
	my @arr=split("\t",$_);
	push @h9hesc2de, $arr[0];
}
close IN;

open IN, "H1DE_PE.editSwitch.txt.4R" or die $!;
while (<IN>) {
	my @arr=split("\t",$_);
	push @h1de2pe, $arr[0];
}
close IN;

open IN, "H9DE_PE.editSwitch.txt.4R" or die $!;
while (<IN>) {
	my @arr=split("\t",$_);
	push @h9de2pe, $arr[0];
}
close IN;

# intersection
my @isect1 = intersect(@h1hesc2de, @h9hesc2de);
my @isect2 = intersect(@h1de2pe, @h9de2pe);

open OUT, ">hESC_to_DE.H1_H9_shared_switching.genelist";
print OUT join ("\n",@isect1), "\n";
close OUT;

open OUT, ">DE_to_PE.H1_H9_shared_switching.genelist";
print OUT join ("\n",@isect2), "\n";
close OUT;







