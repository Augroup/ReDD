#!/bin/bash
# How test_data/rna004 was generated (lab-internal paths; kept for provenance).
#
# Source: HEK293T-WT RNA004 direct RNA (SQK-RNA004, FLO-PRO004RA, run PBG24989), MinKNOW pod5;
#         Dorado 1.1.1 rna004_130bps_sup@v5.2.0 basecalling (--emit-moves), minimap2 -y -ax splice -uf -k14 to hg38,
#         held-out chromosomes chr3/chr11.
#         Test reads = primary alignments in hg38 chr11:500,001-1,000,000 that overlap REDItools candidate A-to-I
#         sites (600 reads, seed 42) whose signal is in the pod5 files available at build time.
# Reference: hg38 chr11:600,001-960,000 renamed "chr11_sub" (coordinates shifted by -600,000).
# NOTE: `pod5 filter` must be run with pod5 <= 0.3.28; files written by newer pod5 versions (read table v6)
#       cannot be opened by Dorado 1.1.1.
set -euo pipefail
POD5_DIR=/scratch/kinfai_root/kinfai0/haorli/20260906_HEK293T_WT_RNA004_pod5
POD5_FILES=$(find "$POD5_DIR" -name "*.pod5" -size +1G | sort)   # completely transferred files only
BAM=/scratch/kinfai_root/kinfai0/haorli/20250815_ReDD_second_revision/20250905_RNA004_data_aligned/aligned/HEK293T-WT/aligned_test_sup.bam
GENOME=/nfs/turbo/umms-kinfai/haorli/20250815_ReDD_second_revision/reference/genome.fa
CAND=/nfs/turbo/umms-kinfai/haorli/20250815_ReDD_second_revision/reference/reditools2_candidates/HEK293T_WT.candidate_sites.tab
NAME=HEK293T-WT_chr11_sub
REGION_START=500001; REGION_END=1000000     # read selection window (hg38 chr11)
OFFSET=600000; REF_END=960000               # mini reference window (hg38 chr11)
OUT=${1:-test_data/rna004}
mkdir -p "$OUT/pod5" "$OUT/reference"

# 1. read ids available in the pod5 files
python - $POD5_FILES <<'EOF'
import sys, pod5
out = open("pod5_have_ids.txt", "w")   # written to the working directory (large)
for fn in sys.argv[1:]:
    with pod5.Reader(fn) as r:
        out.write("".join(f"{x}\n" for x in r.read_ids))
EOF

# 2. reads overlapping candidate sites in the window (signal id = parent read for Dorado split reads, pi tag)
awk -F'\t' -v s=$REGION_START -v e=$REGION_END '$1=="chr11" && $2>=s && $2<=e {print $1"\t"$2-1"\t"$2}' "$CAND" > cand_win.bed
samtools view -L cand_win.bed "$BAM" chr11:$REGION_START-$REGION_END \
  | awk -F'\t' '{id=$1; for(i=12;i<=NF;i++) if($i ~ /^pi:Z:/) id=substr($i,6); print $1"\t"id}' | sort -u > reads_win.txt
awk 'NR==FNR{a[$1];next} ($2 in a)' pod5_have_ids.txt reads_win.txt > reads_sel.txt
cut -f1 reads_sel.txt | sort -u | shuf -n 600 --random-source=<(yes 42) | sort > "$OUT/read_ids.txt"
awk 'NR==FNR{a[$1];next} ($1 in a){print $2}' "$OUT/read_ids.txt" reads_sel.txt | sort -u > signal_ids.txt
samtools view -b -N "$OUT/read_ids.txt" -o test.aligned.bam "$BAM" chr11:$REGION_START-$REGION_END && samtools index test.aligned.bam

# 3. raw signal subset (pod5 <= 0.3.28!)
pod5 filter --ids signal_ids.txt --output "$OUT/pod5/$NAME.pod5" --missing-ok --force-overwrite $POD5_FILES

# 4. Dorado-style unaligned BAM with move tags (alternative input that skips basecalling)
python "$(dirname "$0")/make_test_data.py" test.aligned.bam "$OUT/$NAME.unaligned.bam"

# 5. mini reference and shifted candidate sites
samtools faidx "$GENOME" chr11:$((OFFSET+1))-$REF_END | sed '1s/.*/>chr11_sub/' > "$OUT/reference/chr11_sub.fa"
samtools faidx "$OUT/reference/chr11_sub.fa"
awk -F'\t' -v OFS='\t' -v off=$OFFSET -v end=$REF_END '$1=="chr11" && $2>off && $2<=end {$1="chr11_sub"; $2=$2-off; print}' "$CAND" \
  > "$OUT/$NAME.candidate_sites.tab"
