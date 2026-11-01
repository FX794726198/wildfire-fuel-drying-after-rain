# Wet-fuel days after fire-season rains at ISD stations with measured humidity, compared with ERA5 at the same cell and hour.
# Station T/Td at a fixed synoptic UTC hour near local noon; ERA5 T/Td at the same hour; same ERA5 daily rain (24 h to local noon)
# and the same fire-season months (ERA5 cell climatology) for both, so the only difference is the temperature/humidity source.
import numpy as np, pandas as pd, os, sys
sys.path.insert(0, './code')
from e2_ffmc_global import ffmc, rh
D = '.'; MINCOV = float(os.environ.get('MINCOV', 0.7)); RAIN = os.environ.get('RAIN', 'era5')
S = pd.read_csv(f'{D}/isd/stations_cand.csv', dtype={'USAF': str, 'WBAN': str})
dates = pd.date_range('1979-01-01', '2022-12-31'); mon, yrs = dates.month.values, dates.year.values; nd = len(dates)
E = [np.load(f'{D}/isd/era5s/era5s_{y}.npz') for y in range(1979, 2023)]
ET, ETd = np.concatenate([e['T'] for e in E], 1), np.concatenate([e['Td'] for e in E], 1)
cells = S.cell.values
if RAIN == 'era5':
    P = np.concatenate([np.load(f'{D}/era5g/era5_{y}.npz')['P24'][cells] for y in range(1979, 2023)], 1)
else:                                                    # gauge rain (CPC) on the same cells
    P = np.load(os.environ['RAINNPZ'])['P'][cells].astype(np.float32)
smon = np.load(f'{D}/out_era5/ffmc_obs.npz')['smon'][cells]
def evts(F, Pc, season, thr=10., rt=85.):
    rows = []; ok = np.isfinite(F)
    for i in np.where(season & (Pc >= thr))[0]:
        if i < 3 or i + 31 >= nd or not ok[i - 3:i + 1].all() or not F[i - 1] >= rt: continue
        seg = F[i + 1:i + 31]; kk = np.where(seg >= rt)[0]; rows.append((yrs[i], kk[0] + 1 if len(kk) else 31))
    return np.array(rows).reshape(-1, 2)
# pass 1: QC and gap filling per station; pass 2: FFMC for all stations at once; pass 3: events and statistics
keep, A = [], {k: [] for k in ['T', 'Td', 'Tf', 'Hs', 'eT', 'He']}
for s, r in enumerate(S.itertuples()):
    f = f'{D}/isd/st/{r.USAF}-{r.WBAN}.npz'
    if not os.path.exists(f): continue
    z = np.load(f); T, Td = z['T'].astype(float), z['Td'].astype(float); eT, eTd = ET[s].astype(float), ETd[s].astype(float)
    bad = ~((T > -60) & (T < 60) & (Td > -70) & (Td <= T + 0.5)); T[bad] = np.nan; Td[bad] = np.nan; Td = np.minimum(Td, T)
    # robust QC against ERA5 (differences more than 6 scaled MADs from the station's calendar-month median are dropped)
    for X, Y in ((T, eT), (Td, eTd)):
        for m in range(1, 13):
            k = mon == m; dd = X[k] - Y[k]; med = np.nanmedian(dd); mad = 1.4826 * np.nanmedian(np.abs(dd - med))
            q = np.abs(dd - med) > 6 * max(mad, 0.5); idx = np.where(k)[0][q]; T[idx] = np.nan; Td[idx] = np.nan
    # fill gaps with ERA5 at the same hour plus the station's calendar-month offset
    Tf, Tdf = T.copy(), Td.copy()
    for X, Y in ((Tf, eT), (Tdf, eTd)):
        for m in range(1, 13):
            k = (mon == m); off = np.nanmedian(X[k] - Y[k]); miss = k & ~np.isfinite(X); X[miss] = Y[miss] + (off if np.isfinite(off) else 0)
    keep.append(s)
    for k_, v in zip(['T', 'Td', 'Tf', 'Hs', 'eT', 'He'], [T, Td, Tf, rh(Tf, np.minimum(Tdf, Tf)), eT, rh(eT, eTd)]): A[k_].append(v.astype(np.float32))
A = {k: np.stack(v) for k, v in A.items()}; keep = np.array(keep); Pk = P[keep]; W = np.full(Pk.shape, 12., np.float32)
Fs = ffmc(A['Tf'], A['Hs'], W, Pk); Fe = ffmc(A['eT'], A['He'], W, Pk); print('ffmc done', flush=True)
out = []
for j, s in enumerate(keep):
    r = S.iloc[s]; T, Td, eT, He = A['T'][j].astype(float), A['Td'][j].astype(float), A['eT'][j].astype(float), A['He'][j].astype(float)
    season = smon[s][mon - 1]; have = np.isfinite(T) & np.isfinite(Td)
    cov_e, cov_l = have[season & (yrs < 2001)].mean(), have[season & (yrs >= 2001)].mean()
    es, ee = evts(Fs[j], Pk[j], season), evts(Fe[j], Pk[j], season)
    def dd(e):
        a, b = e[e[:, 0] < 2001, 1], e[e[:, 0] >= 2001, 1]
        return (b.mean() / a.mean() - 1 if len(a) >= 10 and len(b) >= 10 and len(e) >= 30 else np.nan), len(a), len(b)
    ds, ne_s, nl_s = dd(es); de, ne_e, nl_e = dd(ee)
    hm = have & season; se, sl = season & (yrs < 2001), season & (yrs >= 2001)
    Hraw = rh(T, np.minimum(Td, T))
    out.append(dict(sid=f'{r.USAF}-{r.WBAN}', name=r['name'], ctry=r.ctry, lat=r.lat, lon=r.lon, cell=r.cell, hs=int(3 * np.rint((12.5 - r.lon / 15) / 3)),
                    cov_e=cov_e, cov_l=cov_l, d_st=ds, d_era=de, ne_st=ne_s, nl_st=nl_s, ne_era=ne_e, nl_era=nl_e,
                    rec_st_e=es[es[:, 0] < 2001, 1].mean() if len(es) else np.nan, rec_era_e=ee[ee[:, 0] < 2001, 1].mean() if len(ee) else np.nan,
                    # measured-only humidity/temperature change (days with station data in both sources)
                    dRH_st=np.nanmean(Hraw[hm & sl]) - np.nanmean(Hraw[hm & se]), dRH_era=np.nanmean(He[hm & sl]) - np.nanmean(He[hm & se]),
                    dT_st=np.nanmean(T[hm & sl]) - np.nanmean(T[hm & se]), dT_era=np.nanmean(eT[hm & sl]) - np.nanmean(eT[hm & se]),
                    RH_st=np.nanmean(Hraw[hm]), RH_era=np.nanmean(He[hm]),
                    # annual fire-season station-minus-ERA5 RH difference, for break screening
                    **{f'dif{y}': np.nanmean((Hraw - He)[hm & (yrs == y)]) if (hm & (yrs == y)).sum() >= 20 else np.nan for y in range(1979, 2023)}))
R = pd.DataFrame(out); tag = os.environ.get('TAG', 'main'); R.to_csv(f'{D}/res/isd_{tag}.csv', index=False)
k = (R.cov_e >= MINCOV) & (R.cov_l >= MINCOV) & R.d_st.notna() & R.d_era.notna()
print(tag, 'stations', len(R), 'used', k.sum(), 'median d station', round(R.d_st[k].median(), 4), 'ERA5 same hour', round(R.d_era[k].median(), 4),
      'dRH st', round(R.dRH_st[k].median(), 2), 'era', round(R.dRH_era[k].median(), 2))
