#!/usr/bin/perl

use strict;
use warnings;

my @sample=("H1","H9");

foreach my $each (@sample){
my $in=$each."hESC_DE_PE.meanEdit4Isoform.compare";
my $out1=$each."hESC_DE.editSwitch.txt";
my $out2=$each."DE_PE.editSwitch.txt";	
open IN, "$in" or die $!;
open OUT, ">$out1";
open OUT2, ">$out2";
my %gene;
while (<IN>) {
	next if (/^Isoform/);
	my $line=$_;
	my @arr=split("\t",$line);
	push @{$gene{$arr[1]}}, $line;
}

foreach my $key (keys %gene){
	my @each=@{$gene{$key}};
	my $each=@each;
	my ($iso1_s1,$iso1_s2,$iso1_s3,$iso2_s1,$iso2_s2,$iso2_s3);
	if ($each==2) {
		my @arr1=split("\t",$each[0]);
		$iso1_s1=$arr1[3];
		$iso1_s2=$arr1[6];
		$iso1_s3=$arr1[9];

		my @arr2=split("\t",$each[1]);
		$iso2_s1=$arr2[3];
		$iso2_s2=$arr2[6];
		$iso2_s3=$arr2[9];

		next if ($iso1_s1==0 and $iso2_s1==0);
		next if ($iso1_s2==0 and $iso2_s2==0);
		next if ($iso1_s3==0 and $iso2_s3==0);

		my $ratio1=abs(($iso1_s1-$iso1_s2)/$iso1_s1); ## first isoform
		my $ratio2=abs(($iso1_s2-$iso1_s3)/$iso1_s2);

		my $ratio3=abs(($iso2_s1-$iso2_s2)/$iso2_s1); ## second isoform
		my $ratio4=abs(($iso2_s2-$iso2_s3)/$iso2_s2);

		
		if (($iso1_s1<$iso1_s2 and $iso2_s1 >$iso2_s2) or ($iso1_s1>$iso1_s2 and $iso2_s1 <$iso2_s2 )) {
			if ( $ratio1>=0.2 and $ratio2>=0.2 ) {
				print OUT "$arr1[1]\t$arr1[0]\t$arr1[3]\t$arr1[6]\t$arr2[0]\t$arr2[3]\t$arr2[6]\t$ratio1\t$ratio2\n";
			}
		}
		
		if (($iso1_s2<$iso1_s3 and $iso2_s2 >$iso2_s3) or ($iso1_s2>$iso1_s3 and $iso2_s2 <$iso2_s3 )) {
			if ( $ratio3>=0.2 and $ratio4>=0.2 ) {
				print OUT2 "$arr1[1]\t$arr1[0]\t$arr1[6]\t$arr1[9]\t$arr2[0]\t$arr2[6]\t$arr2[9]\t$ratio2\t$ratio3\n";
			}
		}
	}

	elsif($each>2){
		my @sort= sort {(split (/\t/,$b))[3] <=> (split (/\t/,$a))[3] } @each; ## order by the edit level
		my @sort2= sort {(split (/\t/,$b))[6] <=> (split (/\t/,$a))[6] } @each; ## order by the edit level

		my @arr1=split("\t",$sort[0]);
		$iso1_s1=$arr1[3];
		$iso1_s2=$arr1[6];
		$iso1_s3=$arr1[9];

		my @arr2=split("\t",$sort2[0]);
		$iso2_s1=$arr2[3];
		$iso2_s2=$arr2[6];
		$iso2_s3=$arr2[9];

		next if ($arr1[0] eq $arr2[0]); ## if the highest edited isoform remain no change
		next if ($iso1_s1==0 and $iso2_s1==0);
		next if ($iso1_s2==0 and $iso2_s2==0);
		next if ($iso1_s3==0 and $iso2_s3==0);

		my $ratio1=abs(($iso1_s1-$iso1_s2)/$iso1_s1); ## first isoform
		my $ratio2=abs(($iso1_s2-$iso1_s3)/$iso1_s2);

		my $ratio3=abs(($iso2_s1-$iso2_s2)/$iso2_s1); ## second isoform
		my $ratio4=abs(($iso2_s2-$iso2_s3)/$iso2_s2);


		if (($iso1_s1<$iso1_s2 and $iso2_s1 >$iso2_s2) or ($iso1_s1>$iso1_s2 and $iso2_s1 <$iso2_s2 )) {
			if ( $ratio1>=0.2 and $ratio2>=0.2 ) {
				print OUT "$arr1[1]\t$arr1[0]\t$arr1[3]\t$arr1[6]\t$arr2[0]\t$arr2[3]\t$arr2[6]\t$ratio1\t$ratio2\n";
			}
		}
		
		if (($iso1_s2<$iso1_s3 and $iso2_s2 >$iso2_s3) or ($iso1_s2>$iso1_s3 and $iso2_s2 <$iso2_s3 )) {
			if ( $ratio3>=0.2 and $ratio4>=0.2 ) {
				print OUT2 "$arr1[1]\t$arr1[0]\t$arr1[6]\t$arr1[9]\t$arr2[0]\t$arr2[6]\t$arr2[9]\t$ratio2\t$ratio3\n";
			}
		}
	}
	
}

close IN;
close OUT;
close OUT2;


}