#!/bin/bash
# Step 1 (array 0-43): hourly extraction per year. Step 2 (array 0-39): Nelson model per cell chunk (submit after step 1 with --dependency).
#SBATCH --partition=compute-p1
#SBATCH --account=research-ceg-wm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3900M
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/%u/rfnc/logs/nelson_%x_%a.log
B=${RFNC_BASE:-/scratch/$USER/rfnc}; module load 2025 python/3.11.9; source $B/venv/bin/activate; cd $B/code2
if [ "$STEP" = "1" ]; then python n1_hourly.py $((1979+SLURM_ARRAY_TASK_ID)) $B/work/hourly
else python n2_nelson.py $B/work/hourly $B/code2/libdfm.so $SLURM_ARRAY_TASK_ID 40 $B/out/nelson; fi
