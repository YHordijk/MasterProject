#!/bin/bash
#SBATCH -N 2
#SBATCH -t 240:00:00
#SBATCH --ntasks-per-node=64
#SBATCH --partition=tc

echo "== Starting workflow at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"
echo "== Scratch dir. : ${TMPDIR}"

module load ams
module list

$AMSBIN/amspython workflow.py
