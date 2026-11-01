# Forced changes of monthly near-surface temperature and relative humidity from CMIP6 DAMIP, at the ERA5 1.5-deg analysis cells.
# delta_ant = historical(+ssp245 after 2014) - hist-nat ; delta_ghg = hist-GHG - its 1850-1900 mean ; delta_aer = hist-aer - its
# 1850-1900 mean (hist-GHG and hist-aer are single-forcing runs, so their own change is the response; v1.1 wrongly subtracted hist-nat).
# Ensemble mean over members, 11-year centred running mean per calendar month. After the last simulated year: hist-nat held at the mean
# of its last 20 years; hist-GHG and hist-aer extended with the linear trend of their last 20 years. Output npz per model:
# dT, dH [kind, cell, year(1979..2022), month], last simulated year of each experiment.
# Usage: python e6b_damip.py <model> <cells_all.csv> <outdir>
import numpy as np, pandas as pd, xarray as xr, gcsfs, cftime, sys, os
model, cellf, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
d = pd.read_csv(os.environ.get('CAT', 'pangeo.csv')); d = d[(d.table_id == 'Amon') & (d.source_id == model)]
fs = gcsfs.GCSFileSystem(token='anon'); C = pd.read_csv(cellf)
YS = np.arange(1850, 2023); LAST = {}
def series(exp, v, maxmem=10):
    s = d[(d.experiment_id == exp) & (d.variable_id == v)]
    if s.empty: return None
    acc = []
    for mem, g in list(s.groupby('member_id'))[:maxmem]:
        ds = xr.open_zarr(fs.get_mapper(g.sort_values('version').zstore.iloc[-1]), consolidated=True, decode_times=False)
        t = ds.time; tt = cftime.num2date(t.values, t.attrs['units'], t.attrs.get('calendar', 'standard'))
        yr = np.array([x.year for x in tt]); mo = np.array([x.month for x in tt]); k = (yr >= YS[0]) & (yr <= YS[-1])
        la, lo = ds.lat.values, ds.lon.values % 360
        ia = np.abs(la[None, :] - C.lat.values[:, None]).argmin(1)
        io = np.abs(((lo[None, :] - (C.lon.values[:, None] % 360) + 180) % 360) - 180).argmin(1)
        a = ds[v].isel(time=np.where(k)[0]).values[:, ia, io]                   # time x cell
        A = np.full((len(C), len(YS), 12), np.nan, np.float32)
        A[:, yr[k] - YS[0], mo[k] - 1] = a.T
        acc.append(A); LAST[exp] = min(LAST.get(exp, 9999), int(yr[k].max()))
    return np.nanmean(np.stack(acc), 0), len(acc)
def smooth(A):
    B = np.full_like(A, np.nan)
    for j in range(len(YS)):
        w = slice(max(0, j - 5), min(len(YS), j + 6)); B[:, j] = np.nanmean(A[:, w], 1)
    return B
def lastyr(A): return int(YS[np.where(np.isfinite(A).any((0, 2)))[0].max()])
def extend_flat(A):                       # fill years after the last simulated year with the mean of the last 20 years
    L = lastyr(A); k = (YS > L - 20) & (YS <= L); m = np.nanmean(A[:, k], 1); B = A.copy()
    for j in np.where(YS > L)[0]: B[:, j] = np.where(np.isnan(B[:, j]), m, B[:, j])
    return B
def extend_trend(A):                      # fill years after the last simulated year with the linear trend of the last 20 years
    L = lastyr(A); k = (YS > L - 20) & (YS <= L); x = YS[k] - YS[k].mean(); y = A[:, k]
    b = np.nansum((y - np.nanmean(y, 1, keepdims=True)) * x[None, :, None], 1) / (x ** 2).sum(); m = np.nanmean(y, 1); B = A.copy()
    for j in np.where(YS > L)[0]: B[:, j] = np.where(np.isnan(B[:, j]), m + b * (YS[j] - YS[k].mean()), B[:, j])
    return B
res = {}; info = {}
for v in ['tas', 'hurs']:
    h, nh = series('historical', v); s, ns = series('ssp245', v)
    h = np.where(np.isnan(h), s, h)                                        # 2015-2022 from ssp245
    n, nn = series('hist-nat', v); n = extend_flat(n)
    kinds = {'ant': h - n}
    g = series('hist-GHG', v); a = series('hist-aer', v)
    base = (YS >= 1850) & (YS <= 1900)
    if g is not None: G = extend_trend(g[0]); kinds['ghg'] = G - np.nanmean(G[:, base], 1, keepdims=True)
    if a is not None: A_ = extend_trend(a[0]); kinds['aer'] = A_ - np.nanmean(A_[:, base], 1, keepdims=True)
    for k_, x in kinds.items(): res[(k_, v)] = smooth(x)[:, YS >= 1979]
    info[v] = dict(hist=nh, ssp=ns, nat=nn, ghg=g[1] if g else 0, aer=a[1] if a else 0)
kinds = sorted({k for k, _ in res if (k, 'tas') in res and (k, 'hurs') in res})
np.savez_compressed(f'{out}/damip_{model}.npz', kinds=np.array(kinds), dT=np.stack([res[(k, 'tas')] for k in kinds]),
                    dH=np.stack([res[(k, 'hurs')] for k in kinds]), years=YS[YS >= 1979], last=str(LAST))
print(model, kinds, info, LAST, {k: (float(np.nanmean(res[(k, 'tas')][:, -11:]) - np.nanmean(res[(k, 'tas')][:, :11])),
                                float(np.nanmean(res[(k, 'hurs')][:, -11:]) - np.nanmean(res[(k, 'hurs')][:, :11]))) for k in kinds})
