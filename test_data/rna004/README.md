# RNA004 test data

HEK293T-WT RNA004 direct-RNA reads (SQK-RNA004 / FLO-PRO004RA, MinKNOW pod5, Dorado `sup` basecalling) that overlap REDItools A-to-I candidate
sites in hg38 chr11:500,001-1,000,000.

| file | content |
|---|---|
| `pod5/HEK293T-WT_chr11_sub.pod5` | MinKNOW pod5 subset (`pod5 filter`) of the test reads (Dorado split reads are stored under their parent read id, as in a real run) |
| `HEK293T-WT_chr11_sub.unaligned.bam` | Dorado-style unaligned BAM with move tags (`mv`, `ts`, `ns`, `pi`, `sp`) — optional input that skips basecalling |
| `reference/chr11_sub.fa` | hg38 chr11:600,001-960,000 as contig `chr11_sub` (position = hg38 position − 600,000) |
| `HEK293T-WT_chr11_sub.candidate_sites.tab` | HEK293T-WT candidate sites in that window, coordinates shifted accordingly |
| `read_ids.txt` | the test read ids |
| `make_test_data.sh`, `make_test_data.py` | how the files were produced (lab-internal paths) |

Run the pipeline on it (environment `ReDD_RNA004` activated, model weights in `scripts/models/rna004/`):
```
bash test_data/rna004/run_test.sh [output_dir] [CPU|GPU] [pod5|bam]
```
