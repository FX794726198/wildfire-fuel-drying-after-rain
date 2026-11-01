# Batch jobs on the DelftBlue cluster: CMIP6 runs, fire detections and rain products

Base folder on the cluster: `/scratch/$USER/rfnc` (set `RFNC_BASE` to change it), with `code/`, `out/`, `work/`, `logs/` and `venv/`.

1. Copy this folder, together with the current `code/*.py`, to `$RFNC_BASE/code` and run `bash setup.sh` on the login node (needs internet).
2. Submit the jobs:
   - `sbatch job_cmip6.sh`: 57 tasks (24 CMIP6 models x historical+SSP2-4.5 / historical+SSP5-8.5 / piControl), reading the Pangeo CMIP6 Zarr store on Google Cloud; about 10-60 min each.
   - `sbatch job_fire.sh`: FIRMS MODIS and VIIRS yearly all-country archives (https://firms.modaps.eosdis.nasa.gov/data/country/) and GFED4.1s daily emissions.
   - `sbatch job_rain.sh`: CPC global gauge analysis (0.5°, 1979-2022) and CHIRPS (0.25°, 1981-2022).
3. When the jobs have finished, run `bash collect.sh` and copy `rfnc_results_cmip6.tgz`, `fire_*.npz`, `rain_*.npz` and `logs_tail.txt` to the local working folder.
