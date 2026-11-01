# Run the global ERA5 wet-fuel-day analysis for several input variants.
# Usage: python e2_run.py <era5 dir> <out dir> <variant> [rain npz]
#   variants: obs | warm (RH trend removed) | dry (T trend removed) | neither | windobs | rt90 | txrn (Tmax + RHmin, CMIP-comparable)
#             | rain:<name> (obs T/RH with daily rain from <rain npz>, key P, same cells x days)
import numpy as np, pandas as pd, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2_ffmc_global import ffmc, rh, detrend, events, summarise, norain_rec
src, out, var = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
Y0, Y1, SPLIT = int(os.environ.get('Y0', 1979)), int(os.environ.get('Y1', 2022)), int(os.environ.get('SPLIT', 2001))
t0 = time.time()
Z = [np.load(f'{src}/era5_{y}.npz') for y in range(Y0, Y1 + 1)]
cat = lambda k: np.concatenate([z[k] for z in Z], 1)
dates = pd.date_range(f'{Y0}-01-01', f'{Y1}-12-31'); mon, yrs = dates.month.values, dates.year.values
P = cat('P24')
if var == 'txrn': T, H = cat('Tmax'), cat('RHmin')
else: T = cat('T12'); H = rh(T, cat('Td12'))
W = cat('W12') if var == 'windobs' else np.full(T.shape, 12., np.float32)
kind = var.split(':')[2] if var.startswith('rain:') and var.count(':') == 2 else var
if kind in ('warm', 'neither'): H = np.clip(detrend(H, mon, yrs, 1990), 1, 100)
if kind in ('dry', 'neither'): T = detrend(T, mon, yrs, 1990)
if var.startswith('rain:'):
    P = np.load(sys.argv[4])['P'].astype(np.float32); assert P.shape == T.shape
if var.startswith('delta:'):                     # remove a DAMIP forced change: delta:<model>:<kind>[:T|:H]
    _, mdl, kind = var.split(':')[:3]; which = var.split(':')[3] if var.count(':') == 3 else 'TH'
    z = np.load(f"{os.environ.get('DAMIP', './damip')}/damip_{mdl}.npz"); ki = list(z['kinds']).index(kind)
    yi = yrs - int(z['years'][0])
    if 'T' in which: T = T - np.nan_to_num(z['dT'][ki][:, yi, mon - 1])
    if 'H' in which: H = np.clip(H - np.nan_to_num(z['dH'][ki][:, yi, mon - 1]), 1, 100)
rt = 90. if var == 'rt90' else 85.
F = ffmc(T, H, W, P); print('ffmc done', round(time.time() - t0))
ev, smon = events(F, P, mon, yrs, rt=rt)
if var == 'obs':
    ev['rec_nr'] = norain_rec(ev, F, T, H, W)
    np.savez_compressed(f'{out}/ffmc_obs.npz', F=F.astype(np.float16), smon=smon)
elif os.environ.get('SAVEF'):
    np.savez_compressed(f"{out}/ffmc_{var.replace(':', '_')}.npz", F=F.astype(np.float16), smon=smon)
R = summarise(ev, F, P, H, T, smon, mon, yrs, SPLIT, Y0, Y1)
tag = var.replace(':', '_')
ev.to_csv(f'{out}/events_{tag}.csv.gz', index=False); R.to_csv(f'{out}/cells_{tag}.csv', index=False)
print(tag, 'cells', len(R), 'events', len(ev), 'median d', round(R.d.median(), 4), 'shorter', round((R.d < 0).mean(), 3),
      'dRH', round((R.H_l - R.H_e).median(), 2), 'dT', round((R.T_l - R.T_e).median(), 2), 'sec', round(time.time() - t0))
