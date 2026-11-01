# Drying time after fire-season rains with daily maximum temperature and minimum relative humidity (the inputs that CMIP6 provides)
# measured at ISD stations, compared with ERA5 sampled at the same 3-hourly synoptic hours of the same local days.
# Two station versions:
#  A (ERA5-assisted, as in v1.1): QC against ERA5, gaps filled with ERA5 plus the station's monthly offset; each source has its own
#    events (pre-rain FFMC >= 85 in that source).
#  B (station-only): QC against the station's own monthly median (6 scaled MADs), no gap filling; events need valid station FFMC on
#    the 7 days before the rain and every day until the fuels are dry again; ERA5 is evaluated on exactly the same events.
# Same ERA5 daily rain (24 h to local noon) and the same observed fire-season months (ERA5 noon climatology) for all versions.
# Drying times capped at 15 days and pooled over events, as for the CMIP6 comparison.  Output res/isd_txrn.csv
# Fire-season T/RH changes: dRH_st/dT_st from the version-A measured (ERA5-screened) days; dRH_stB/dT_stB from the version-B values
# (station-only screening, unfilled) on all valid fire-season days, each with ERA5 on the same days.
import numpy as np, pandas as pd, os, sys
sys.path.insert(0, './code')
from e2_ffmc_global import ffmc
D = '.'; CAP = 15
S = pd.read_csv(f'{D}/isd/stations_cand.csv', dtype={'USAF': str, 'WBAN': str})
dates = pd.date_range('1979-01-01', '2022-12-31'); mon, yrs = dates.month.values, dates.year.values; nd = len(dates)
E = [np.load(f'{D}/isd/era5t/era5t_{y}.npz') for y in range(1979, 2023)]
ETx, EHn = np.concatenate([e['Tx'] for e in E], 1), np.concatenate([e['Hn'] for e in E], 1)
cells = S.cell.values
P = np.concatenate([np.load(f'{D}/era5g/era5_{y}.npz')['P24'][cells] for y in range(1979, 2023)], 1)
smon = np.load(f'{D}/out_era5/ffmc_obs.npz')['smon'][cells]
def mad_qc(X, ref=None):
    """drop values > 6 scaled MADs from the calendar-month median of X (or of X - ref)."""
    X = X.copy()
    for m in range(1, 13):
        k = mon == m; dd = X[k] - (ref[k] if ref is not None else 0); med = np.nanmedian(dd)
        mad = 1.4826 * np.nanmedian(np.abs(dd - med)); bad = np.abs(dd - med) > 6 * max(mad, 0.5)
        idx = np.where(k)[0][bad]; X[idx] = np.nan
    return X
A = {k: [] for k in ['TA', 'HA', 'TB', 'HB', 'ET', 'EH', 'Tm', 'Hm']}; keep = []
for s, r in enumerate(S.itertuples()):
    f = f'{D}/isd/st3/{r.USAF}-{r.WBAN}.npz'
    if not os.path.exists(f): continue
    z = np.load(f); ok = z['n3'] >= 7
    Tx = np.where(ok, z['Tx'], np.nan).astype(float); Hn = np.where(ok, z['Hn'], np.nan).astype(float)
    eT, eH = ETx[s].astype(float), EHn[s].astype(float)
    # A: QC against ERA5 and ERA5-based gap filling
    TA, HA = mad_qc(Tx, eT), mad_qc(Hn, eH); both = np.isfinite(TA) & np.isfinite(HA); TA[~both] = np.nan; HA[~both] = np.nan
    Tm, Hm = TA.copy(), HA.copy()                                    # measured (QC-passed) values for the humidity trend
    for X, Y in ((TA, eT), (HA, eH)):
        for m in range(1, 13):
            k = mon == m; off = np.nanmedian(X[k] - Y[k]); miss = k & ~np.isfinite(X); X[miss] = Y[miss] + (off if np.isfinite(off) else 0)
    # B: QC against the station's own climatology, no filling
    TB, HB = mad_qc(Tx), mad_qc(Hn); both = np.isfinite(TB) & np.isfinite(HB); TB[~both] = np.nan; HB[~both] = np.nan
    keep.append(s)
    for k_, v in zip(A, [TA, np.clip(HA, 1, 100), TB, HB, eT, eH, Tm, Hm]): A[k_].append(v.astype(np.float32))
A = {k: np.stack(v) for k, v in A.items()}; keep = np.array(keep); Pk = P[keep]; W = np.full(Pk.shape, 12., np.float32)
FA = ffmc(A['TA'], A['HA'], W, Pk); FB = ffmc(A['TB'], A['HB'], W, Pk); FE = ffmc(A['ET'], A['EH'], W, Pk); print('ffmc done', flush=True)
def evts(F, Pc, season, pre=3, need=None):
    rows = []; ok = np.isfinite(F) if need is None else need
    for i in np.where(season & (Pc >= 10.))[0]:
        if i < pre or i + 31 >= nd or not ok[i - pre:i + 1].all() or not F[i - 1] >= 85: continue
        seg = F[i + 1:i + 31]; kk = np.where(seg >= 85)[0]; rec = kk[0] + 1 if len(kk) else 31
        if not ok[i + 1:i + 1 + min(rec, 30)].all(): continue
        rows.append((i, rec))
    return np.array(rows, int).reshape(-1, 2)
def rec_at(F, idx):
    out = []
    for i in idx:
        seg = F[i + 1:i + 31]; kk = np.where(seg >= 85)[0]; out.append(kk[0] + 1 if len(kk) else 31)
    return np.array(out, int)
def change(i, rec):
    if len(i) == 0: return np.nan, 0, 0
    late = yrs[i] >= 2001; r = np.minimum(rec, CAP); a, b = r[~late], r[late]
    return (b.mean() / a.mean() - 1 if len(a) >= 10 and len(b) >= 10 and len(r) >= 30 else np.nan), len(a), len(b)
def trendB(j, season):
    # version B: fire-season changes of daily max T and min RH from the station-only (own-climatology QC, unfilled) values on all
    # valid fire-season days, and ERA5 at the same hours on exactly the same days
    m = np.isfinite(A['TB'][j]) & np.isfinite(A['HB'][j]) & season; se, sl = m & (yrs < 2001), m & (yrs >= 2001)
    ne, nl = (season & (yrs < 2001)).sum(), (season & (yrs >= 2001)).sum()
    f = lambda X: np.nanmean(X[sl]) - np.nanmean(X[se]) if se.any() and sl.any() else np.nan
    return dict(covB_e=se.sum() / ne if ne else 0, covB_l=sl.sum() / nl if nl else 0, dRH_stB=f(A['HB'][j]), dRH_eraB=f(A['EH'][j]),
                dT_stB=f(A['TB'][j]), dT_eraB=f(A['ET'][j]))
out = []
for j, s in enumerate(keep):
    r = S.iloc[s]; season = smon[s][mon - 1]
    eA, eE = evts(FA[j], Pk[j], season), evts(FE[j], Pk[j], season)
    eB = evts(FB[j], Pk[j], season, pre=7, need=np.isfinite(A['TB'][j]) & np.isfinite(A['HB'][j]))
    dA, neA, nlA = change(eA[:, 0], eA[:, 1]); dE, neE, nlE = change(eE[:, 0], eE[:, 1])
    dB, neB, nlB = change(eB[:, 0], eB[:, 1]); dBE, _, _ = change(eB[:, 0], rec_at(FE[j], eB[:, 0]))
    m = np.isfinite(A['Tm'][j]) & season; se, sl = m & (yrs < 2001), m & (yrs >= 2001)
    cov_e, cov_l = m[season & (yrs < 2001)].mean() if (season & (yrs < 2001)).any() else 0, m[season & (yrs >= 2001)].mean() if (season & (yrs >= 2001)).any() else 0
    out.append(dict(sid=f'{r.USAF}-{r.WBAN}', ctry=r.ctry, lat=r.lat, lon=r.lon, cell=r.cell, cov_e=cov_e, cov_l=cov_l,
                    dA_st=dA, dA_era=dE, neA=neA, nlA=nlA, dB_st=dB, dB_era=dBE, neB=neB, nlB=nlB,
                    dRH_st=np.nanmean(A['Hm'][j][sl]) - np.nanmean(A['Hm'][j][se]), dRH_era=np.nanmean(A['EH'][j][sl & np.isfinite(A['Hm'][j])]) - np.nanmean(A['EH'][j][se & np.isfinite(A['Hm'][j])]),
                    dT_st=np.nanmean(A['Tm'][j][sl]) - np.nanmean(A['Tm'][j][se]), dT_era=np.nanmean(A['ET'][j][sl & np.isfinite(A['Tm'][j])]) - np.nanmean(A['ET'][j][se & np.isfinite(A['Tm'][j])]),
                    **trendB(j, season)))
R = pd.DataFrame(out); R.to_csv(f'{D}/res/isd_txrn.csv', index=False)
okA = R.dA_st.notna() & R.dA_era.notna() & (R.cov_e >= 0.7) & (R.cov_l >= 0.7); okB = R.dB_st.notna() & R.dB_era.notna()
print('A stations', okA.sum(), 'median st', round(R.dA_st[okA].median() * 100, 2), 'era', round(R.dA_era[okA].median() * 100, 2),
      'dRH st', round(R.dRH_st[okA].median(), 2), 'era', round(R.dRH_era[okA].median(), 2))
print('B stations', okB.sum(), 'median st', round(R.dB_st[okB].median() * 100, 2), 'era same events', round(R.dB_era[okB].median() * 100, 2),
      'events', int(R.neB[okB].sum() + R.nlB[okB].sum()), 'dRH stB', round(R.dRH_stB[okB].median(), 2), 'eraB', round(R.dRH_eraB[okB].median(), 2),
      'dT stB', round(R.dT_stB[okB].median(), 2), 'eraB', round(R.dT_eraB[okB].median(), 2))
