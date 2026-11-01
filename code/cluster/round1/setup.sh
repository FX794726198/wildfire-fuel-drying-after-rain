#!/bin/bash
# one-time: venv for the NC-upgrade jobs
set -e
B=${RFNC_BASE:-/scratch/$USER/rfnc}; mkdir -p $B/out $B/logs $B/work
module load 2025 python/3.11.9 2>/dev/null || module load 2024r1 python
python3 -m venv $B/venv && source $B/venv/bin/activate
pip install -q --upgrade pip
pip install -q numpy pandas scipy "xarray" "zarr<3" gcsfs cftime numcodecs requests h5py netCDF4
python -c "import xarray, zarr, gcsfs, cftime, h5py, netCDF4; print('venv ok', xarray.__version__, zarr.__version__)"
