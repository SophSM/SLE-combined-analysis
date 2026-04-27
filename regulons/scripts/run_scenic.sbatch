#!/bin/bash
#SBATCH -J scenic_sle_100runs
#SBATCH -p defq
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 48:00:00
#SBATCH -o /mnt/data/amedina/ssalazar/regulons_SLE/scripts/logs/head_%j.out
#SBATCH -e /mnt/data/amedina/ssalazar/regulons_SLE/scripts/logs/head_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sogusama02@gmail.com

set -euo pipefail

BASE="/mnt/data/amedina/ssalazar/regulons_SLE"
SCRIPTS="${BASE}/scripts"
LOGS="${SCRIPTS}/logs"
BIN="${BASE}/bin"

mkdir -p "${LOGS}" "${BASE}/work" "${BASE}/out" "${BASE}/singularity" "${BASE}/tmp" "${BASE}/.nextflow_home"
cd "${BASE}"

# ---------- Deterministic runtime ----------
export PATH="${BIN}:${PATH}"

export NXF_HOME="${BASE}/.nextflow_home"
export NXF_ANSI_LOG=false
export NXF_PROGRESS=false
export NXF_OPTS="-Xms512m -Xmx6g -Djava.io.tmpdir=${BASE}/tmp"
export NXF_SINGULARITY_CACHEDIR="${BASE}/singularity"

# Nextflow 21.04.3: log file via env var (NOT -log)
export NXF_LOG_FILE="${LOGS}/nextflow_${SLURM_JOB_ID}.log"

# ---------- Modules ----------
if [ -f /etc/profile.d/modules.sh ]; then
  . /etc/profile.d/modules.sh
fi

module purge || true
module load singularity/3.7.0 || module load singularity || true
module load oracle-java/8u371

# ---------- Sanity checks ----------
echo "[$(date)] Host: $(hostname)"
echo "[$(date)] PWD:  $(pwd)"
echo "[$(date)] java -version:"
java -version 2>&1
echo "[$(date)] nextflow path: ${BIN}/nextflow"
"${BIN}/nextflow" -version

# ---------- Run ----------
"${BIN}/nextflow" \
  -C "${SCRIPTS}/sle_scenic.config" \
  run vib-singlecell-nf/vsn-pipelines \
  -entry scenic \
  -r v0.25.0 \
  -resume \
  -work-dir "${BASE}/work" \
  -qs 1
