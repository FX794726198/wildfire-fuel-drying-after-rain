# Data for the supplementary figures S1-S4 (metric illustration, robustness maps, trend-removal maps, fire rate by FFMC).
import numpy as np, pandas as pd, sys, os, glob
B = './'; O = B + 'out_era5/'; F = B + 'figs/data/'
C = pd.read_csv(B + 'cells_all.csv'); SEL = set(pd.read_csv(B + 'cells_sel.csv').cell)
def cells(tag):
    R = pd.read_csv(O + f'cells_{tag}.csv'); R = R[R.cell.isin(SEL)].merge(C[['cell', 'lat', 'lon']], on='cell'); R['rel'] = 100 * R.d; return R
# ---- S1: composite FFMC around rain events (early vs late), distribution of drying time, example cell
Z = np.load(O + 'ffmc_obs.npz'); FF = Z['F'].astype(np.float32); smon = Z['smon']
E = pd.read_csv(O + 'events_obs.csv.gz'); E = E[E.cell.isin(SEL)].copy(); E['late'] = (E.y >= 2001).astype(int)
lags = np.arange(-3, 11); nd = FF.shape[1]
E = E[(E.i + 10 < nd)]
comp = FF[E.cell.values[:, None], E.i.values[:, None] + lags[None, :]]
K = pd.DataFrame(comp, columns=lags); K['cell'] = E.cell.values; K['late'] = E.late.values
cm = K.groupby(['cell', 'late']).mean()                                   # mean per cell and period, then over cells
rows = []
rng = np.random.default_rng(0); cl = np.array(sorted(SEL))
for late in (0, 1):
    m = cm.xs(late, level='late'); bs = [m.loc[m.index.intersection(rng.choice(cl, len(cl)))].mean().values for _ in range(300)]
    lo, hi = np.nanpercentile(bs, 2.5, 0), np.nanpercentile(bs, 97.5, 0)
    for j, l in enumerate(lags): rows.append(dict(late=late, lag=l, ffmc=m[l].mean(), lo=lo[j], hi=hi[j]))
pd.DataFrame(rows).to_csv(F + 'figS1_composite.csv', index=False)
h = E.assign(r=E.rec.clip(upper=16)).groupby(['late', 'r']).size().unstack(0).fillna(0); h = h / h.sum()
h.to_csv(F + 'figS1_hist.csv')
# example cell: one in the western United States with many events; one fire season
ex = E[(E.cell.isin(C[(C.lat.between(36, 42)) & (C.lon.between(-123, -118))].cell))].groupby('cell').size().idxmax()
yrs = E[E.cell == ex].groupby('y').size(); yex = int(yrs.idxmax())
d = pd.date_range('1979-01-01', '2022-12-31'); k = np.where(d.year == yex)[0]
P = np.load(B + f'era5g/era5_{yex}.npz')['P24'][ex]
pd.DataFrame(dict(date=d[k], ffmc=FF[ex, k], rain=P, season=smon[ex, d.month[k] - 1])).to_csv(F + 'figS1_example.csv', index=False)
E[(E.cell == ex) & (E.y == yex)][['i', 'rec', 'p']].assign(date=lambda q: d[q.i]).to_csv(F + 'figS1_example_events.csv', index=False)
print('example cell', ex, C.loc[ex, ['lat', 'lon']].values, 'year', yex)
# ---- S2: maps for alternative inputs
M = cells('obs')[['cell', 'lat', 'lon', 'rel']].rename(columns={'rel': 'era5'})
for tag, nm in [('rain_chirps', 'chirps'), ('txrn', 'txrn'), ('windobs', 'wind')]:
    M = M.merge(cells(tag)[['cell', 'rel']].rename(columns={'rel': nm}), on='cell', how='left')
NE = pd.read_csv(B + 'res/nelson_events.csv.gz')
for v, t, nm in [('m1_14', 0.10, 'nelson1'), ('m10_14', 0.12, 'nelson10')]:
    g = NE[(NE['var'] == v) & (NE.thr == t)].assign(per=lambda q: np.where(q.y < 2001, 'e', 'l')).groupby(['cell', 'per']).agg(r=('rc', 'mean'), n=('rc', 'size')).unstack()
    g = g[(g[('n', 'e')] >= 10) & (g[('n', 'l')] >= 10)]; s = 100 * (g[('r', 'l')] / g[('r', 'e')] - 1)
    M = M.merge(s.rename(nm).reset_index(), on='cell', how='left')
M.to_csv(F + 'figS2_maps.csv', index=False); print(M.drop(columns=['cell', 'lat', 'lon']).median().round(2).to_dict())
# ---- S3: trend-removal maps and temperature change
T = cells('obs')[['cell', 'lat', 'lon']].assign(dT=lambda q: 0)
o = cells('obs'); T['dT'] = (o.T_l - o.T_e).values
for tag in ['warm', 'dry', 'neither']: T = T.merge(cells(tag)[['cell', 'rel']].rename(columns={'rel': tag}), on='cell', how='left')
T.to_csv(F + 'figS3_maps.csv', index=False); print(T[['dT', 'warm', 'dry', 'neither']].median().round(2).to_dict())
# ---- S4: relative fire rate by FFMC class (within cell-year-month strata) and by days since rain
Z5 = np.load(B + 'dbo/fire_modis.npz'); N = Z5['N']; off = (pd.Timestamp('2001-01-01') - d[0]).days; nf = nd - off
Nf = N[:, :nf]; Fw = FF[:, off:]; mon = d.month.values[off:]; yr = d.year.values[off:]
A = pd.read_csv(B + 'res/fire_MH_modis_cells.csv'); cl = A.cell.values
bins = [0, 60, 70, 75, 80, 83, 85, 87, 89, 91, 93, 95, 101]
num = np.zeros(len(bins) - 1); den = np.zeros(len(bins) - 1); per_reg = []
for c in cl:
    s = smon[c, mon - 1] & np.isfinite(Fw[c])
    df = pd.DataFrame(dict(st=yr[s] * 100 + mon[s], f=Fw[c, s], n=Nf[c, s]))
    df['exp'] = df.groupby('st').n.transform('mean'); df['b'] = np.clip(np.digitize(df.f, bins) - 1, 0, len(bins) - 2)
    g = df.groupby('b')[['n', 'exp']].sum(); num[g.index] += g.n.values; den[g.index] += g.exp.values
pd.DataFrame(dict(lo=bins[:-1], hi=bins[1:], ratio=num / den, fires=num)).to_csv(F + 'figS4_ffmc_bins.csv', index=False)
print(pd.read_csv(F + 'figS4_ffmc_bins.csv').round(3).to_string())
