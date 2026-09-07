# ReDD RNA004 pipeline image: latest code from the GitHub repository, conda environment (incl. the vendored
# uncalled4 ReDD fork), Dorado + RNA004 basecalling model, ReDD RNA004 model weights and the bundled test data.
#
#   docker build -t redd-rna004 --build-arg CACHEBUST=$(date +%s) .           # CACHEBUST forces a fresh git clone
#   docker build -t redd-rna004 --build-arg REDD_BRANCH=<branch or tag> .     # another revision
#
#   docker run --rm --gpus all -v /path/to/data:/data redd-rna004 \
#       bash test_data/rna004/run_test.sh /data/redd_test GPU                  # smoke test (CPU: drop --gpus all, use CPU)
#   docker run --rm --gpus all -v /path/to/data:/data redd-rna004 bash -c \
#       "python generate_script.py rna004 --pipeline_mode bash --input_pod5 /data/pod5 --ref_genome /data/genome.fa \
#           --output_path /data/run1 --output_name sample1 --device GPU && bash /data/run1/run.pbs"
#
# Torch and Dorado ship their own CUDA libraries; the host only needs the NVIDIA driver and nvidia-container-toolkit.
FROM condaforge/miniforge3:24.11.3-2

ARG REDD_REPO=https://github.com/Augroup/ReDD.git
ARG REDD_BRANCH=main
ARG CACHEBUST=1

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends git curl ca-certificates procps \
 && rm -rf /var/lib/apt/lists/*

# 1. latest pipeline code (CACHEBUST only invalidates the layer cache)
RUN echo "cachebust=${CACHEBUST}" \
 && git clone --depth 1 --branch "${REDD_BRANCH}" "${REDD_REPO}" /opt/ReDD
WORKDIR /opt/ReDD

# 2. conda environment (pip section installs ./uncalled4_ReDD)
RUN mamba env create -n ReDD_RNA004 --file environment_rna004.yaml \
 && mamba clean -afy \
 && /opt/conda/envs/ReDD_RNA004/bin/pip cache purge \
 && rm -rf uncalled4_ReDD/build

# 3. Dorado + RNA004 basecalling model (-> /opt/ReDD/software/dorado), 4. model weights (-> scripts/models/rna004)
RUN bash scripts/rna004/install_dorado.sh \
 && bash scripts/rna004/download_model.sh

ENV PATH=/opt/conda/envs/ReDD_RNA004/bin:/opt/ReDD/software/dorado/bin:$PATH \
    CONDA_DEFAULT_ENV=ReDD_RNA004 \
    CONDA_PREFIX=/opt/conda/envs/ReDD_RNA004
RUN echo 'source /opt/conda/etc/profile.d/conda.sh && conda activate ReDD_RNA004' > /etc/profile.d/redd.sh

CMD ["bash"]
