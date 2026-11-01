# Drying times after a FIXED set of observed rain events under modified weather (paired comparison; review of v1.1).
# For each variant the FFMC is recomputed from the modified daily weather, and the drying time is taken after exactly the same rain
# events (cell, day) as in the reference run, so that differences between variants come only from the weather, not from a change in
# which rains qualify or which cells are analysed. Output fx/rec_<ref>_<variant>.npy (int8 drying time per reference event).
# Usage: python e7_fixed.py <ref> <variant> [<variant> ...]
#   ref: noon (ERA5 noon T/RH, ERA5 rain) | chirps (noon T/RH, CHIRPS rain) | txrn (ERA5 Tmax/RHmin, ERA5 rain, observed noon fire season)
#   variants: obs | warm | dry | neither | dmp:<model>:<kind>[:T|:H] (DAMIP forced change removed, damip2/) |
#             hsc:<f>[:far] (RH trend scaled by f; ':far' only in cells more than 3 deg from an ISD station)
import numpy as np, pandas as pd, sys, os, time
sys.path.insert(0, './code')
from e2_ffmc_global import ffmc, rh, detrend, events
D = '.'; ref, VAR = sys.argv[1], sys.argv[2:]
SEL = pd.read_csv(f'{D}/code/cells_sel.csv').cell.values                  # 1,704 analysis cells (global ids)
Y0 = 1981 if ref == 'chirps' else 1979                                     # CHIRPS starts in 1981
dates = pd.date_range(f'{Y0}-01-01', '2022-12-31'); mon, yrs = dates.month.values, dates.year.values
Z = [np.load(f'{D}/era5g/era5_{y}.npz') for y in range(Y0, 2023)]
cat = lambda k: np.concatenate([z[k][SEL] for z in Z], 1)
P = cat('P24') if ref != 'chirps' else np.load(f'{D}/dbo/rain_chirps.npz')['P'][SEL].astype(np.float32)
if ref == 'txrn': T0, H0 = cat('Tmax'), cat('RHmin')
else: T0 = cat('T12'); H0 = rh(T0, cat('Td12'))
del Z
W = np.full(T0.shape, 12., np.float32); row = {c: j for j, c in enumerate(SEL)}
smon_obs = np.load(f'{D}/out_era5/ffmc_obs.npz')['smon'][SEL]
evf = f'{D}/fx/events_{ref}.csv.gz'
if not os.path.exists(evf):                                                  # reference event set
    if ref == 'noon': ev = pd.read_csv(f'{D}/out_era5/events_obs.csv.gz')
    elif ref == 'chirps': ev = pd.read_csv(f'{D}/out_era5/events_rain_chirps.csv.gz')
    else:
        F = ffmc(T0, H0, W, P); ev, _ = events(F, P, mon, yrs, smon=smon_obs); ev['cell'] = SEL[ev.cell.values]
    ev = ev[ev.cell.isin(set(SEL))].reset_index(drop=True); ev[['cell', 'y', 'm', 'i', 'rec', 'p']].to_csv(evf, index=False)
ev = pd.read_csv(evf); ci, ii = ev.cell.map(row).values, ev.i.values
far = None
def recs(F):
    seg = F[ci[:, None], ii[:, None] + np.arange(1, 31)[None, :]]               # event x 30 following days
    hit = seg >= 85.; first = hit.argmax(1) + 1; first[~hit.any(1)] = 31; return first.astype(np.int8)
for var in VAR:
    out = f'{D}/fx/rec_{ref}_{var.replace(":", "_")}.npy'
    if os.path.exists(out): continue
    t0 = time.time(); T, H = T0, H0
    if var in ('warm', 'neither'): H = np.clip(detrend(H0, mon, yrs, 1990), 1, 100)
    if var in ('dry', 'neither'): T = detrend(T0, mon, yrs, 1990)
    if var.startswith('dmp:'):
        _, mdl, kind = var.split(':')[:3]; which = var.split(':')[3] if var.count(':') == 3 else 'TH'
        z = np.load(f'{D}/damip2/damip_{mdl}.npz'); ki = list(z['kinds']).index(kind); yi = yrs - int(z['years'][0])   # noon/txrn refs start in 1979
        if 'T' in which: T = T0 - np.nan_to_num(z['dT'][ki][SEL][:, yi, mon - 1])
        if 'H' in which: H = np.clip(H0 - np.nan_to_num(z['dH'][ki][SEL][:, yi, mon - 1]), 1, 100)
    if var.startswith('hsc:'):
        f = float(var.split(':')[1]); Hd = detrend(H0, mon, yrs, 1990); H = np.clip(f * H0 + (1 - f) * Hd, 1, 100)
        if var.endswith(':far'):
            if far is None:
                S = pd.read_csv(f'{D}/res/isd_main.csv'); S = S[S.d_st.notna()]; C = pd.read_csv(f'{D}/code/cells_sel.csv')
                dist = np.sqrt((C.lat.values[:, None] - S.lat.values[None, :]) ** 2 +
                               (((C.lon.values[:, None] - S.lon.values[None, :] + 180) % 360 - 180) * np.cos(np.radians(C.lat.values[:, None]))) ** 2)
                far = dist.min(1) > 3
            H = np.where(far[:, None], H, H0)
    F = ffmc(T, H, W, P); r = recs(F); np.save(out, r)
    print(ref, var, 'events', len(r), 'mean rec', round(float(r.mean()), 3), 'sec', round(time.time() - t0), flush=True)
