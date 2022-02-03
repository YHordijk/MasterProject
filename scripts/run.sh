#!/bin/bash
#SBATCH -N 1
#SBATCH -t 240:00:00
#SBATCH --ntasks 16
#SBATCH --partition=tc

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"
echo "== Scratch dir. : ${TMPDIR}"

module load ams
module list

$AMSBIN/amspython job_runner.py > $SLURM_JOBID.log
