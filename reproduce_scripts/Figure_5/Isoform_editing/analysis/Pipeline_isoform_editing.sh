## This pipeline is used to identify A-to-I editing on isoform level
## Author: Bo Li 05-10-2024

# First step: filter out the low-confidence editing sites identified by REDD

perl 01_filterEditSite.pl

# Second step: calculate editing levels for each isoform

perl 02_callMeanEditLevelIsoform.pl

# third step: calculate editing diversity

perl 03_editDiversityEachSample.pl

# fourth step: calculate editing diversity at site level

perl 04_EditDiversityWithSite.pl

# fifth step: analysis editing switching 

perl 05_editSwitch.pl

# sixth step: compare and summerize the editing switching results

perl 06_compare_H1_H9_switching_genes.pl