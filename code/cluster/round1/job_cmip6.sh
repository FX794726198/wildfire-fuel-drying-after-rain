#!/bin/bash
#SBATCH --job-name=rf_cmip6
#SBATCH --partition=compute-p1
#SBATCH --account=research-ceg-wm
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=56G
#SBATCH --array=1-57%20
#SBATCH --output=/scratch/%u/rfnc/logs/cmip6_%a.log
B=${RFNC_BASE:-/scratch/$USER/rfnc}; module load 2025 python/3.11.9 2>/dev/null || module load 2024r1 python; source $B/venv/bin/activate
cd $B/code
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" cmip6_jobs.csv); M=${LINE%,*}; S=${LINE#*,}
CAT=$B/work/pangeo-cmip6.csv python e3_cmip6.py $M $S cells_all.csv $B/out/cmip6
