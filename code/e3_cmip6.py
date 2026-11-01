# CMIP6 daily (Pangeo GCS zarr): FFMC from tasmax + hursmin (standard noon proxy) + pr, fixed 12 km/h wind, at the ERA5 1.5-deg
# analysis cells (nearest model grid point). Output: per cell x year mean wet-fuel days per rain and event counts.
# Usage: python e3_cmip6.py <model> <scenario> <cells.csv> <outdir>
#   scenario: hist245 | hist585 | histnat | picontrol      (member: first available r*i1p1f*, or env MEMBER)
import numpy as np, pandas as pd, xarray as xr, gcsfs, sys, os, time, cftime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2_ffmc_global import ffmc, events
model, scen, cellf, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
CAT = os.environ.get('CAT', 'pangeo-cmip6.csv')
if not os.path.exists(CAT): os.system(f'curl -s -o {CAT} https://storage.googleapis.com/cmip6/pangeo-cmip6.csv')
d = pd.read_csv(CAT); d = d[(d.table_id == 'day') & (d.source_id == model)]
fs = gcsfs.GCSFileSystem(token='anon')
HV = os.environ.get('HV', 'hursmin')
EXP = {'hist245': ['historical', 'ssp245'], 'hist585': ['historical', 'ssp585'], 'histnat': ['hist-nat'], 'picontrol': ['piControl']}[scen]
VARS = ['tasmax', HV, 'pr']
def members(exp):
    s = d[d.experiment_id == exp]; ok = [m for m, g in s.groupby('member_id') if set(VARS) <= set(g.variable_id)]
    return sorted(ok, key=lambda m: (('i1p1' not in m), int(m.split('i')[0][1:])))
mem = os.environ.get('MEMBER') or [m for m in members(EXP[0]) if all(m in members(e) for e in EXP)][0]
C = pd.read_csv(cellf)
def load(exp, v):
    z = d[(d.experiment_id == exp) & (d.member_id == mem) & (d.variable_id == v)].sort_values('version').zstore.iloc[-1]
    ds = xr.open_zarr(fs.get_mapper(z), consolidated=True, decode_times=False)
    t = ds.time; tt = cftime.num2date(t.values, t.attrs['units'], t.attrs.get('calendar', 'standard'))
    la, lo = ds.lat.values, ds.lon.values % 360
    ia = np.abs(la[None, :] - C.lat.values[:, None]).argmin(1)
    io = np.abs(((lo[None, :] - (C.lon.values[:, None] % 360) + 180) % 360) - 180).argmin(1)
    yr = np.array([x.year for x in tt]); mo = np.array([x.month for x in tt]); dy = np.array([x.day for x in tt])
    if scen == 'hist245' or scen == 'hist585':
        k = (yr >= 1950) if exp == 'historical' else (yr <= 2100)
    elif scen == 'picontrol': k = yr < yr.min() + int(os.environ.get('NPI', 250))
    else: k = yr >= 1950
    if os.environ.get('YMAX'): k = k & (yr <= int(os.environ['YMAX']))      # large-ensemble runs stop at 2022
    kk = np.where(k)[0]; out = np.empty((len(C), len(kk)), np.float32); B = int(os.environ.get('BLK', 730))
    for j in range(0, len(kk), B):                                # read time blocks, keep only the analysis cells
        a = ds[v].isel(time=slice(kk[j], kk[min(j + B, len(kk)) - 1] + 1)).values
        out[:, j:j + a.shape[0]] = a[:, ia, io].T
    return out, yr[k], mo[k], dy[k]
t0 = time.time(); X = {v: [] for v in VARS}; YR, MO = [], []
for e in EXP:
    L = {v: load(e, v) for v in VARS}                              # align the variables on common days (some stores miss days)
    keys = [L[v][1] * 10000 + L[v][2] * 100 + L[v][3] for v in VARS]; com = keys[0]
    for kk_ in keys[1:]: com = np.intersect1d(com, kk_)
    for v, kk_ in zip(VARS, keys):
        _, ix = np.unique(kk_, return_index=True); pos = ix[np.searchsorted(kk_[ix], com)]; X[v].append(L[v][0][:, pos])
    j0 = keys[0]; _, ix = np.unique(j0, return_index=True); pos = ix[np.searchsorted(j0[ix], com)]
    YR.append(L[VARS[0]][1][pos]); MO.append(L[VARS[0]][2][pos]); print(model, e, mem, 'loaded', round(time.time() - t0), 'days', len(com), flush=True)
T = np.concatenate(X['tasmax'], 1) - 273.15; H = np.clip(np.concatenate(X[HV], 1), 1, 100); P = np.concatenate(X['pr'], 1) * 86400.
yr, mo = np.concatenate(YR), np.concatenate(MO)
F = ffmc(T, H, np.full(T.shape, 12., np.float32), P)
ref = (1950, 2014) if scen != 'picontrol' else (int(yr.min()), int(yr.max()))
SM = np.load(os.environ['SMON'])['smon'] if os.environ.get('SMON') else None   # fixed (observed) fire-season months [cell, 12]
ev, smon = events(F, P, mo, yr, ref=ref, smon=SM)
g = ev.groupby(['cell', 'y']).agg(rec=('rec', 'mean'), rec15=('rec', lambda x: x.clip(upper=15).mean()), n=('rec', 'size'), p=('p', 'mean')).reset_index()
# fire-season mean noon proxies per cell-year (for scaling with humidity)
seas = smon[:, mo - 1]
yy = np.unique(yr); HT = []
for y in yy:
    k = yr == y; m = seas[:, k]
    HT.append(pd.DataFrame(dict(cell=np.arange(len(C)), y=y, H=np.nanmean(np.where(m, H[:, k], np.nan), 1), T=np.nanmean(np.where(m, T[:, k], np.nan), 1),
                                wet=(np.where(m, F[:, k], np.nan) < 85).sum(1), rd=(np.where(m, P[:, k], 0) >= 1).sum(1))))
HT = pd.concat(HT); g = HT.merge(g, on=['cell', 'y'], how='left')
g.to_csv(f'{out}/{model}_{scen}_{mem}.csv.gz', index=False)
if os.environ.get('SMON'): np.save(f'{out}/{model}_{scen}_{mem}_mF.npy', np.stack([np.nanmean(F[:, mo == m], 1) for m in range(1, 13)], 1).astype(np.float32))
print(model, scen, mem, 'cells with events', ev.cell.nunique(), 'events', len(ev), 'sec', round(time.time() - t0))
