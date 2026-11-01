#!/bin/bash
#SBATCH --job-name=rf_fire
#SBATCH --partition=compute-p1
#SBATCH --account=research-ceg-wm
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --array=0-2
#SBATCH --output=/scratch/%u/rfnc/logs/fire_%a.log
B=${RFNC_BASE:-/scratch/$USER/rfnc}; module load 2025 python/3.11.9 2>/dev/null || module load 2024r1 python; source $B/venv/bin/activate
cd $B/code; SRC=(modis viirs gfed)
python e4_fire.py ${SRC[$SLURM_ARRAY_TASK_ID]} cells_all.csv $B/work/fire
