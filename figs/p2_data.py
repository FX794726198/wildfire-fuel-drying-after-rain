# Data for Fig. 2 (fires and wet fuels): Mantel-Haenszel rate ratios per cell and region, dose-response with cell bootstrap.
import numpy as np, pandas as pd, sys, os, json
B = './'; F = B + 'figs/data/'; tag = sys.argv[1] if len(sys.argv) > 1 else 'modis'
C = pd.read_csv(B + 'cells_all.csv')
A = pd.read_csv(B + f'res/fire_MH_{tag}_cells.csv').merge(C[['cell', 'lat', 'lon']], on='cell'); A.to_csv(F + f'fig2_cells_{tag}.csv', index=False)
rng = np.random.default_rng(3); rows = []
for r, g in A.groupby('ar6'):
    if len(g) < 15: continue
    bs = []
    for _ in range(1000):
        k = g.sample(len(g), replace=True, random_state=int(rng.integers(1e9))); bs.append(k.num.sum() / k.den.sum())
    rows.append(dict(ar6=r, n=len(g), rr=g.num.sum() / g.den.sum(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5)))
pd.DataFrame(rows).sort_values('rr').to_csv(F + f'fig2_regions_{tag}.csv', index=False)
E = pd.read_csv(B + f'res/fire_events_{tag}.csv.gz'); cells = E.cell.unique(); Gc = {c: g for c, g in E.groupby('cell')}
def curve(e): s = e.groupby('rc')[['fev', 'fctl']].sum(); return s.fev / s.fctl
base = curve(E); bs = []
for _ in range(500):
    bs.append(curve(pd.concat([Gc[c] for c in rng.choice(cells, len(cells))])).reindex(base.index))
bs = np.array(bs)
D = pd.DataFrame(dict(rc=base.index, ratio=base.values, lo=np.nanpercentile(bs, 2.5, 0), hi=np.nanpercentile(bs, 97.5, 0), n=E.groupby('rc').size().values))
D.to_csv(F + f'fig2_dose_{tag}.csv', index=False); print(D.round(3).to_string())
S = json.load(open(B + f'res/fire_summary_{tag}.json')); print(S)
