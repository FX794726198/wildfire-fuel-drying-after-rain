# Daily rain on the ERA5 1.5-deg analysis cells from gauge-based products (rain-source robustness of the global analysis).
#  cpc    : NOAA CPC Unified Gauge-Based Global Daily Precipitation 0.5 deg, 1979-2022 (mean of the 3x3 0.5-deg boxes inside each 1.5-deg cell)
#  chirps : CHIRPS v2.0 daily p25 (50S-50N), 1981-2022 (mean of the 0.25-deg boxes inside each cell; cells outside -> NaN)
# Output npz key P [cell, day], same day axis as e1/e2 for the given years. Usage: python e5_rain.py <cpc|chirps> <cells_all.csv> <workdir>
import numpy as np, pandas as pd, xarray as xr, sys, os, subprocess
src, cellf, wd = sys.argv[1:4]; os.makedirs(wd, exist_ok=True)
C = pd.read_csv(cellf)
y0 = 1979 if src == 'cpc' else 1981; Y1 = 2022
dates = pd.date_range(f'{y0}-01-01', f'{Y1}-12-31'); P = np.full((len(C), len(dates)), np.nan, np.float32)
def wget(url, f):
    if not os.path.exists(f): subprocess.run(['wget', '-q', '-c', '-O', f + '.part', url], check=True); os.rename(f + '.part', f)
for y in range(y0, Y1 + 1):
    if src == 'cpc':
        f = f'{wd}/precip.{y}.nc'; wget(f'https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/precip.{y}.nc', f); v = 'precip'
    else:
        f = f'{wd}/chirps-v2.0.{y}.days_p25.nc'
        if not os.path.exists(f) and os.path.exists(f'{os.environ.get("CHIRPS_RAW", "chirps_raw")}/chirps-v2.0.{y}.days_p25.nc'): f = f'{os.environ.get("CHIRPS_RAW", "chirps_raw")}/chirps-v2.0.{y}.days_p25.nc'
        else: wget(f'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p25/chirps-v2.0.{y}.days_p25.nc', f)
        v = 'precip'
    ds = xr.open_dataset(f); lat = ds.lat.values; lon = ds.lon.values % 360
    a = ds[v].values.astype(np.float32); a[a < 0] = np.nan                     # time x lat x lon
    out = np.full((len(C), a.shape[0]), np.nan, np.float32)
    for c, la, lo in zip(C.cell, C.lat, C.lon):
        ka = np.abs(lat - la) < 0.75; ko = np.abs(((lon - lo % 360 + 180) % 360) - 180) < 0.75
        if ka.any() and ko.any():
            x = a[:, ka][:, :, ko].reshape(a.shape[0], -1)
            good = np.isfinite(x).mean(0) > 0.5
            if good.any(): out[c] = np.nanmean(x[:, good], 1)
    i = (pd.Timestamp(y, 1, 1) - dates[0]).days; P[:, i:i + out.shape[1]] = out
    print(src, y, np.nanmean(out), flush=True)
np.savez_compressed(f'{wd}/rain_{src}.npz', P=P, y0=y0)
