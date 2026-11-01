# Global daily fire activity on the ERA5 1.5-deg analysis cells.
#  modis|viirs : NASA FIRMS yearly all-country archives (MODIS C6.1 2001-2024, VIIRS S-NPP 375 m 2012-2024); presumed vegetation fires
#                (type 0), MODIS confidence >= 30 / VIIRS confidence n|h; counts and FRP summed per cell and LOCAL solar day.
#  gfed        : GFED4.1s monthly burned area and C emissions x daily emission fractions (2003-2023, 2017+ beta), summed to 1.5 deg.
# Usage: python e4_fire.py <modis|viirs|gfed> <cells_all.csv> <workdir>
import numpy as np, pandas as pd, sys, os, zipfile, io, glob, subprocess
src, cellf, wd = sys.argv[1:4]; os.makedirs(wd, exist_ok=True)
C = pd.read_csv(cellf); key = {(int(round((la + 90) / 1.5)), int(round((lo % 360) / 1.5)) % 240): c for c, la, lo in zip(C.cell, C.lat, C.lon)}
lut = np.full((121, 240), -1, np.int32)
for (a, b), c in key.items(): lut[a, b] = c
def cellof(lat, lon): return lut[np.rint((lat + 90) / 1.5).astype(int).clip(0, 120), (np.rint((lon % 360) / 1.5).astype(int)) % 240]
def wget(url, f):
    if not os.path.exists(f): subprocess.run(['wget', '-q', '-c', '-O', f + '.part', url], check=True); os.rename(f + '.part', f)
if src in ('modis', 'viirs'):
    y0 = 2001 if src == 'modis' else 2012; dates = pd.date_range(f'{y0}-01-01', '2024-12-31'); d0 = dates[0]
    N = np.zeros((len(C), len(dates)), np.float32); FRP = np.zeros_like(N)
    for y in range(y0, 2025):
        name = f'modis_{y}_all_countries.zip' if src == 'modis' else f'viirs-snpp_{y}_all_countries.zip'
        f = f'{wd}/{name}'; wget(f'https://firms.modaps.eosdis.nasa.gov/data/country/zips/{name}', f)
        with zipfile.ZipFile(f) as z:
            for nm in z.namelist():
                if not nm.endswith('.csv'): continue
                for d in pd.read_csv(z.open(nm), usecols=['latitude', 'longitude', 'acq_date', 'acq_time', 'confidence', 'type', 'frp'], chunksize=2_000_000):
                    d = d[d.type == 0]
                    d = d[d.confidence.astype(str).isin(['n', 'h', 'nominal', 'high'])] if src == 'viirs' else d[pd.to_numeric(d.confidence, errors='coerce') >= 30]
                    t = pd.to_datetime(d.acq_date) + pd.to_timedelta(d.acq_time // 100 * 60 + d.acq_time % 100, unit='m') + pd.to_timedelta(d.longitude / 15, unit='h')
                    di = (t.dt.floor('D') - d0).dt.days.values; c = cellof(d.latitude.values, d.longitude.values)
                    k = (c >= 0) & (di >= 0) & (di < len(dates))
                    np.add.at(N, (c[k], di[k]), 1); np.add.at(FRP, (c[k], di[k]), d.frp.values[k])
        print(src, y, N.sum(), flush=True)
    np.savez_compressed(f'{wd}/fire_{src}.npz', N=N, FRP=FRP, y0=y0)
else:
    import h5py
    dates = pd.date_range('2003-01-01', '2023-12-31'); BA = np.zeros((len(C), len(dates)), np.float32); EC = np.zeros_like(BA)
    la = np.arange(89.875, -90, -0.25); lo = np.arange(-179.875, 180, 0.25)
    cc = cellof(np.repeat(la, 1440), np.tile(lo, 720)).reshape(720, 1440)
    for y in range(2003, 2024):
        name = f'GFED4.1s_{y}.hdf5' if y <= 2016 else f'GFED4.1s_{y}_beta.hdf5'
        f = f'{wd}/{name}'; wget(f'https://www.geo.vu.nl/~gwerf/GFED/GFED4/{name}', f)
        with h5py.File(f) as h:
            area = h['ancill/grid_cell_area'][:]
            for m in range(1, 13):
                g = h[f'emissions/{m:02d}']; ba = h[f'burned_area/{m:02d}/burned_fraction'][:] * area / 1e6 if f'burned_area/{m:02d}' in h else np.zeros_like(area)
                ec = g['C'][:] * area                                   # g C per month
                nd = pd.Period(f'{y}-{m:02d}').days_in_month
                for dd in range(1, nd + 1):
                    fr = g['daily_fraction'][f'day_{dd}'][:]
                    i = (pd.Timestamp(y, m, dd) - dates[0]).days; k = cc >= 0
                    BA[:, i] = np.bincount(cc[k], (ba * fr)[k], len(C)); EC[:, i] = np.bincount(cc[k], (ec * fr)[k], len(C))
        print('gfed', y, BA.sum(), flush=True)
    np.savez_compressed(f'{wd}/fire_gfed.npz', BA=BA, EC=EC, y0=2003)
