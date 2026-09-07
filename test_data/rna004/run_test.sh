#!/bin/bash
# Smoke test of the RNA004 pipeline, starting from pod5 (Dorado basecalling included), on a small set of
# HEK293T-WT reads (hg38 chr11:600,001-960,000 subset).
#
# Usage (from an activated ReDD_RNA004 environment, Dorado installed via scripts/rna004/install_dorado.sh):
#   bash test_data/rna004/run_test.sh [output_dir] [CPU|GPU] [pod5|bam]
#     pod5 (default): basecall the test pod5 with Dorado, then run the pipeline
#     bam           : skip basecalling and use the provided Dorado-style BAM (with move tags)
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
DATA=$REPO/test_data/rna004
OUT=${1:-$REPO/test_data/rna004_results}
DEVICE=${2:-CPU}
MODE=${3:-pod5}
NAME=HEK293T-WT_chr11_sub

READS_ARG=()
[ "$MODE" = "bam" ] && READS_ARG=(--input_bam "$DATA/$NAME.unaligned.bam")

python "$REPO/generate_script.py" rna004 \
    --pipeline_mode bash \
    --output_name "$NAME" \
    --output_path "$OUT" \
    --input_pod5 "$DATA/pod5" \
    "${READS_ARG[@]}" \
    --ref_genome "$DATA/reference/chr11_sub.fa" \
    --ref_candidate_sites "$DATA/$NAME.candidate_sites.tab" \
    --device "$DEVICE" \
    --num_split 1 --threads_basecall 4 --threads_extract 4 --threads_predict 4 --max_cores_resources 8 \
    --coverage_cutoff 5 --ratio_cutoff 0.1 \
    --conda_env "${CONDA_PREFIX:-ReDD_RNA004}"

cd "$OUT"
bash run.pbs
echo "--- snakemake log tail ---"
tail -n 5 "REDD_logs/REDD_$NAME.log"
echo "--- outputs ---"
ls -la outputs
echo "molecule-level predictions: $(wc -l < outputs/$NAME.prediction.genome.txt)"
echo "sites (coverage>=5):        $(wc -l < outputs/$NAME.site.bed)"
echo "filtered sites:             $(wc -l < outputs/$NAME.flt.genome.tab)"
