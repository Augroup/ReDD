#!/usr/bin/perl
use strict;
use warnings;


## calculate the editing diversity of gene isoform

open SUM, ">EditingDiversity.sharedEditSite.sum";
print SUM "cell_line\ttotalSite\tcount_site_as_same\tcount_site_as_diff\tcount_gene_as_same\tcount_gene_as_dif\n";

my $cutoff=0.1; ## ratio for determine editing 

my @files=glob("*.ALUAnno");
foreach my $each (@files){
	$each=~/(.*).ALUAnno/;
	my $sample=$1;
	my %gene;
	open IN, "$each" or die $!;
	open OUT, ">$each.I_NI.txt";
	while (<IN>) {
		chomp;
		my @arr=split("\t",$_);
		my $pos=$arr[6]."\t".$arr[7]."\t".$arr[8];  ## genomic position with gene ID
		push @{$gene{$pos}}, $arr[4]; 
	}

	foreach my $key (keys %gene){
	    my @edit=@{$gene{$key}};
	    my $edit=@edit;
	    my $i=0;
	    my $ni=0;
	    if ($edit>=2) { # this genomic position is shared by at least two isoforms
	    	foreach my $a (@edit){
	    		if ($a>=$cutoff) {
	    			$i++;
	    		}
	    		else{
	    			$ni++;
	    		}
	    	}

	    	print OUT "$key\t$i\t$ni\n";
	    }
	}

	close IN;
	close OUT;

	open IN, "$each.I_NI.txt" or die $!;
	my $totalsite=0;
	my $countsitesame=0;
	my $countsitedif=0;
	#my $countgenesame=0;
	#my $countgenedif=0;
	my @genesame=();
	my @genedif=();
	while (<IN>) {
		chomp;
		$totalsite++;
		my @arr=split("\t",$_);
		if ($arr[3]>0 and $arr[4]>0) {
			$countsitedif++;
			push @genedif,$arr[2];
		}
		else{
			$countsitesame++;
			push @genesame,$arr[2];
		}
	}

	my @genesame_u=uniq(@genesame);
	my @genedif_u=uniq(@genedif);

	my $genesame_u=@genesame_u;
	my $genedif_u=@genedif_u;

	print SUM "$sample\t$totalsite\t$countsitesame\t$countsitedif\t$genesame_u\t$genedif_u\n";

}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
