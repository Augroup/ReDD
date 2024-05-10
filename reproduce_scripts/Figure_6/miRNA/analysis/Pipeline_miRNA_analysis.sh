## This pipeline is used to identify the impact of A-to-I editing on miRNA-mRNA interaction
## Author: Bo Li 05-10-2024

# First step: identify significantly differentially expressed gene isoforms and get editing level for these isoforms

perl 01_EditLevel4DEG.pl

# Second step: run miRanDa to predict the miRNA binding sites for each isoform before and after A-to-I editing

perl 02_run_miRanDa.pl

# Third step: summarize the results

perl 03_summerize_miRanDa.pl