# Installation
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
mamba create -n AIediting --file environment.yaml
```
# Usage
```
source activate AIediting

cd workflow
snakemake -c 6 -p cache/AFG-H1_ONT_directRNA_1.hdf5
```
