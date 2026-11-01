# CMIP6 summary step on the DelftBlue cluster

1. Copy this folder to `/scratch/$USER/rfnc/code3`.
2. After all CMIP6 tasks of `round1/job_cmip6.sh` have finished, on a compute node (one task, 8 GB, 2 h):
   `cd /scratch/$USER/rfnc/code3 && source ../venv/bin/activate && python a5_cmip6.py /scratch/$USER/rfnc/out/cmip6 cmip6`
   This writes `res/cmip6_hist.csv`, `cmip6_proj.csv`, `cmip6_pi.csv`, `cmip6_series.csv` and `cmip6_cells.csv.gz`.
3. Copy `code3/res/` to the local working folder, together with `fire_viirs.npz` and `rain_cpc.npz`.
4. Split any file larger than 380 MB before copying (`split -b 300m`, with the md5 of the whole file in a `.md5` file).
