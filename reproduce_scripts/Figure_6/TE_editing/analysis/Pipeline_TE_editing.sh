## This pipeline is used to analyze the A-to-I editing on transposable elements
## Author: Bo Li 05-10-2024

## Step one: prepare editing level input files

perl 00_callMeanEditLevel_samplespecific.pl

## Step two: compare the A-to-I editing levels between TE and non-TE transcripts

perl 01_TE_vs_nonTE_cellline.pl

## Step three: analyze the A-to-I editing among TE type and TE family

perl 02_TEtype.cmp.pl

perl 02_TEfamily.cmp.pl

## Step four: analyze the A-to-I editing for TE-gene chimeric transcripts

perl 03_TE_Gene_edit.pl

## Step five: analyze the A-to-I editing among different forms of ALU

perl 04_compareEditLevelDifAluForm.pl