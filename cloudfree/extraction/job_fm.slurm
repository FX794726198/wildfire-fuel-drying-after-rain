#!/bin/bash
# Example SLURM job for extract_fm.py on 16 cores (2 satellites x 26 tiles x ~920 granules). Set W to your working folder,
# which holds code/ (this folder, with cells_cloudfree.csv and cells_all.csv) and an empty out/.
#SBATCH --job-name=firemask
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3800M
W=${W:-$PWD}
export TMPDIR=${TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}; mkdir -p $TMPDIR
python $W/code/extract_fm.py $W/code/urls.csv $W/code/cells_cloudfree.csv $W/code/cells_all.csv $W/out 16
