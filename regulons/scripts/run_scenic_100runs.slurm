#!/bin/bash
#SBATCH -J scenic_sle_100runs
#SBATCH -p defq
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 3-00:00:00
#SBATCH -o /mnt/data/amedina/dramirez/Reg_network/final_scenic/scripts/logs/head_%j.out
#SBATCH -e /mnt/data/amedina/dramirez/Reg_network/final_scenic/scripts/logs/head_%j.err
#SBATCH --no-requeue

set -euo pipefail

BASE="/mnt/data/amedina/dramirez/Reg_network/final_scenic"
SCRIPTS="${BASE}/scripts"
LOGS="${SCRIPTS}/logs"
mkdir -p "${LOGS}" "${BASE}/work" "${BASE}/out" "${BASE}/singularity"

echo "[$(date)] Host: $(hostname)"
echo "[$(date)] PWD:  ${BASE}"

module purge
module load oracle-java/8u371
module load singularity/3.7.0

export PATH="${BASE}/bin:${PATH}"
export NXF_VER=21.04.3

echo "[$(date)] java -version:"
java -version || true
echo "[$(date)] nextflow path: $(which nextflow)"
nextflow -version

cd "${BASE}"

# Run ID para logs (mismo jobid ya lo tenemos con %j, esto ayuda cuando corres interactivo)
RUNSTAMP="$(date +%Y%m%d_%H%M%S)"

# Ejecuta Nextflow y guarda log
nextflow -C "${SCRIPTS}/sle_scenic.config" \
  run vib-singlecell-nf/vsn-pipelines \
  -entry scenic \
  -r v0.25.0 \
  -resume \
  -work-dir "${BASE}/work" \
  -with-trace "${LOGS}/execution_trace.txt" \
  -with-report "${LOGS}/execution_report.html" \
  -with-timeline "${LOGS}/execution_timeline.html" \
  -with-dag "${LOGS}/pipeline_dag.svg" \
  2>&1 | tee "${LOGS}/nextflow_${SLURM_JOB_ID}_${RUNSTAMP}.log"

echo "[$(date)] Done."
