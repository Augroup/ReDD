#!/bin/bash
# Download the RNA004 ReDD model weights (not stored in git) into scripts/models/rna004/.
#   bash scripts/rna004/download_model.sh [model_name]     (default: general)
set -euo pipefail
MODEL=${1:-general}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
DEST=$REPO/scripts/models/rna004
mkdir -p "$DEST"
URL=https://reddexamples.s3.us-east-2.amazonaws.com/RNA004/$MODEL.pt
echo "Downloading $URL"
curl -L -o "$DEST/$MODEL.pt" "$URL"
if [ -f "$DEST/$MODEL.pt.sha256" ]; then
    (cd "$DEST" && sha256sum -c "$MODEL.pt.sha256")
fi
echo "Saved to $DEST/$MODEL.pt"
