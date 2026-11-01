# Like-for-like humidity and temperature comparison for the DAMIP part (review of v1.1): fire-season change 1979-2000 -> 2001-2022 of
# ERA5 daily-mean (all-hour) temperature and relative humidity vs the DAMIP anthropogenic response in monthly-mean (all-hour) fields,
# in the 1,704 analysis cells and the observed fire-season months; also ERA5 noon values.  Output res/damip_rh.csv
import numpy as np, pandas as pd, glob, os
D = '.'; SEL = pd.read_csv(f'{D}/code/cells_sel.csv').cell.values
smon = np.load(f'{D}/out_era5/ffmc_obs.npz')['smon'][SEL]                       # cells x 12
dates = pd.date_range('1979-01-01', '2022-12-31'); mon, yrs = dates.month.values, dates.year.values
season = smon[:, mon - 1]; early, late = yrs < 2001, yrs >= 2001
def dchg(X):
    a = np.where(season & early, X, np.nan); b = np.where(season & late, X, np.nan)
    return np.nanmean(b, 1) - np.nanmean(a, 1)
Z = [np.load(f'{D}/era5m/era5m_{y}.npz') for y in range(1979, 2023)]
Tm, Hm = np.concatenate([z['Tm'] for z in Z], 1), np.concatenate([z['Hm'] for z in Z], 1)
G = [np.load(f'{D}/era5g/era5_{y}.npz') for y in range(1979, 2023)]
import sys; sys.path.insert(0, f'{D}/code'); from e2_ffmc_global import rh
T12 = np.concatenate([g['T12'][SEL] for g in G], 1); H12 = rh(T12, np.concatenate([g['Td12'][SEL] for g in G], 1))
rows = [dict(source='ERA5 daily mean', dT=np.median(dchg(Tm)), dRH=np.median(dchg(Hm))),
        dict(source='ERA5 noon', dT=np.median(dchg(T12)), dRH=np.median(dchg(H12)))]
# monthly fire-season weights: each (cell, year, month) weighted by its number of days, as for the daily means
yy = np.arange(1979, 2023); mm = np.arange(1, 13); ndays = np.array([[pd.Period(f'{y}-{m:02d}').days_in_month for m in mm] for y in yy])
for f in sorted(glob.glob(f'{D}/damip2/damip_*.npz')):
    z = np.load(f); mdl = os.path.basename(f)[6:-4]; Y = z['years']; iy = np.searchsorted(Y, yy)
    for k, kind in enumerate(z['kinds']):
        for v, key in (('dT', 'dT'), ('dRH', 'dH')):
            X = z[key][k][SEL][:, iy, :]                                               # cells x years x 12
            w = smon[:, None, :] * ndays[None]; e, l = yy < 2001, yy >= 2001
            ch = (np.nansum(X[:, l] * w[:, l], (1, 2)) / w[:, l].sum((1, 2))) - (np.nansum(X[:, e] * w[:, e], (1, 2)) / w[:, e].sum((1, 2)))
            rows.append(dict(source=mdl, kind=kind, var=v, val=np.nanmedian(ch)))
R = pd.DataFrame(rows); R.to_csv(f'{D}/res/damip_rh.csv', index=False)
P = R[R.kind.notna()].pivot_table(index=['source'], columns=['kind', 'var'], values='val'); print(R[R.kind.isna()][['source', 'dT', 'dRH']].round(2).to_string()); print(P.round(2).to_string())
