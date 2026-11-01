# Nelson dead fuel moisture model (physical check of the FFMC result)

1. Copy this folder to `/scratch/$USER/rfnc/code2` and run `bash build_nelson.sh` (clones NFDRS4 from https://github.com/firelab/NFDRS4 and builds `libdfm.so` with g++). Then `python test_nelson.py code2/libdfm.so` (10-h moisture should jump after day 30 and dry back over several days).
2. Submit the two steps:
   `J1=$(sbatch --parsable --job-name=n1 --array=0-43%15 --export=ALL,STEP=1 job_nelson.sh)`
   `sbatch --job-name=n2 --array=0-39 --dependency=afterok:$J1 --export=ALL,STEP=2 job_nelson.sh`
   (step 2 holds about 43 cells x 385,704 hours x 4 variables in memory; raise `--cpus-per-task` to 8 if it runs out of memory).
3. Copy `/scratch/$USER/rfnc/out/nelson/*.npz` to the local working folder and run `code/a6_nelson.py`.
