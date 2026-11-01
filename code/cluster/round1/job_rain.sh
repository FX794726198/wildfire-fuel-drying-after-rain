#!/bin/bash
#SBATCH --job-name=rf_rain
#SBATCH --partition=compute-p1
#SBATCH --account=research-ceg-wm
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --array=0-1
#SBATCH --output=/scratch/%u/rfnc/logs/rain_%a.log
B=${RFNC_BASE:-/scratch/$USER/rfnc}; module load 2025 python/3.11.9 2>/dev/null || module load 2024r1 python; source $B/venv/bin/activate
cd $B/code; SRC=(cpc chirps)
python e5_rain.py ${SRC[$SLURM_ARRAY_TASK_ID]} cells_all.csv $B/work/rain
