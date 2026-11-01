# 95% cell-bootstrap intervals for the four fire checks of Fig. 2d (Mantel-Haenszel ratio, median over cells,
# 10 days after rain vs matched dry spells, gradient-boosting observed/predicted ratio refitted as in a4_fire.py).
# Usage (from the nc/ folder): python code/a4b_ci.py dbo/fire_modis.npz modis
import numpy as np, pandas as pd, sys, json
fire_f, tag = sys.argv[1], sys.argv[2]; NB = 1000; rng = np.random.default_rng(7); out = {}
A = pd.read_csv(f'res/fire_MH_{tag}_cells.csv'); n = len(A)
bs = []
for _ in range(NB):
    k = A.iloc[rng.integers(0, n, n)]; bs.append((k.num.sum() / k.den.sum(), np.median(k.num / k.den)))
bs = np.array(bs); out['MH_ci'] = np.percentile(bs[:, 0], [2.5, 97.5]).tolist(); out['MH_cell_median_ci'] = np.percentile(bs[:, 1], [2.5, 97.5]).tolist()
B = pd.read_csv(f'res/fire_events_{tag}.csv.gz'); g = B.groupby('cell')[['fev', 'fctl']].sum(); m = len(g)
bs = [g.iloc[rng.integers(0, m, m)].sum() for _ in range(NB)]; out['avoid_ratio_ci'] = np.percentile([b.fev / b.fctl for b in bs], [2.5, 97.5]).tolist()
# gradient boosting, same design as a4_fire.py, then bootstrap the held-out cells
import lightgbm as lgb
O = 'out_era5/'; Z = np.load(fire_f); N = Z['N']; y0 = int(Z['y0']); C = pd.read_csv('cells_all.csv')
FZ = np.load(O + 'ffmc_obs.npz'); F = FZ['F'].astype(np.float32); smon = FZ['smon']
d79 = pd.date_range('1979-01-01', '2022-12-31'); off = (pd.Timestamp(f'{y0}-01-01') - d79[0]).days
nd = len(d79) - off; Nf = N[:, :nd]; mon = d79.month.values[off:]; yr = d79.year.values[off:]
P = np.concatenate([np.load(f'era5g/era5_{y}.npz')['P24'] for y in range(1979, 2023)], 1)[:, off:]
T12 = np.concatenate([np.load(f'era5g/era5_{y}.npz')['T12'] for y in range(1979, 2023)], 1)[:, off:]
Fw = F[:, off:]
rs = np.random.default_rng(0); sel = rs.choice(A.cell.values, min(600, len(A)), replace=False); X = []
for c in sel:
    s = np.where(smon[c, mon - 1] & np.isfinite(Fw[c]))[0]
    X.append(pd.DataFrame(dict(cell=c, y=yr[s], m=mon[s], T=T12[c, s], P=P[c, s], P1=P[c, np.maximum(s - 1, 0)],
                               lat=C.lat[c], lon=C.lon[c], fy=pd.Series(Nf[c, s]).groupby(yr[s]).transform('mean').values,
                               wet=Fw[c, s] < 85, n=Nf[c, s])))
X = pd.concat(X); feats = ['m', 'T', 'P', 'P1', 'lat', 'lon', 'fy', 'y']; tr = X.y % 2 == 0
mdl = lgb.LGBMRegressor(objective='poisson', n_estimators=400, learning_rate=0.05, num_leaves=63, verbose=-1).fit(X[tr][feats], X[tr].n)
te = X[~tr].assign(pred=mdl.predict(X[~tr][feats]))
h = te.groupby(['cell', 'wet'])[['n', 'pred']].sum().unstack('wet').fillna(0); k = len(h)
r = lambda q: (q[('n', True)].sum() / q[('pred', True)].sum()) / (q[('n', False)].sum() / q[('pred', False)].sum())
out['gbm_ratio'] = float(r(h)); out['gbm_ci'] = np.percentile([r(h.iloc[rng.integers(0, k, k)]) for _ in range(NB)], [2.5, 97.5]).tolist()
json.dump(out, open(f'res/fire_ci_{tag}.json', 'w'), indent=1); print(out)
