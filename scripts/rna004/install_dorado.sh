#!/bin/bash
# Download Dorado and the RNA004 basecalling model used by the ReDD RNA004 pipeline.
#
#   bash scripts/rna004/install_dorado.sh [dorado_version] [install_dir] [model ...]
#
# Defaults: dorado 1.1.1, <repo>/software, model rna004_130bps_sup@v5.2.0
# (the version/model the RNA004 ReDD model was trained and validated with).
# Result: <install_dir>/dorado -> dorado-<version>-linux-x64 (binary at software/dorado/bin/dorado)
#         <install_dir>/dorado/models/<model>
# Dorado is distributed by Oxford Nanopore under its own licence (see the LICENSE file in the tarball);
# it is not bundled with this repository.
set -euo pipefail
VERSION=${1:-1.1.1}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
DEST=${2:-$REPO/software}
shift $(( $# > 2 ? 2 : $# ))
MODELS=("$@")
[ ${#MODELS[@]} -eq 0 ] && MODELS=("rna004_130bps_sup@v5.2.0")

mkdir -p "$DEST"
cd "$DEST"
TARBALL=dorado-${VERSION}-linux-x64.tar.gz
if [ ! -x "dorado-${VERSION}-linux-x64/bin/dorado" ]; then
    echo "Downloading $TARBALL ..."
    curl -L -O "https://cdn.oxfordnanoportal.com/software/analysis/${TARBALL}"
    tar -xzf "$TARBALL"
    rm -f "$TARBALL"
fi
ln -sfn "dorado-${VERSION}-linux-x64" dorado
mkdir -p dorado/models
for m in "${MODELS[@]}"; do
    if [ ! -d "dorado/models/$m" ]; then
        echo "Downloading model $m ..."
        dorado/bin/dorado download --model "$m" --models-directory dorado/models
    fi
done
echo "dorado: $DEST/dorado/bin/dorado ($(dorado/bin/dorado --version 2>&1 | tail -1))"
echo "models: $(ls dorado/models | tr '\n' ' ')"
