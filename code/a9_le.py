# Large-ensemble members: change in drying time (rec15, event-weighted per cell), fire-season dRH and dT between 1979-2000 and
# 2001-2022, median over the analysis cells, as for the 21-model comparison (a5_cmip6.py). Output res/le_members.csv.
import numpy as np, pandas as pd, glob, os
os.chdir('.'); SEL = set(pd.read_csv('cells_sel.csv').cell); CS = pd.read_csv('code/cells_sel.csv').cell.values
def periods(g, a, b):
    x = g[g.n.notna()].assign(s=lambda q: q.rec15 * q.n)
    def pm(r):
        z = x[(x.y >= r[0]) & (x.y <= r[1])]
        k = z.groupby('cell').agg(s=('s', 'sum'), n=('n', 'sum')); k['r'] = k.s / k.n
        h = g[(g.y >= r[0]) & (g.y <= r[1])].groupby('cell')[['H', 'T']].mean()
        return k.join(h)
    j = pm(a).join(pm(b), lsuffix='_e', rsuffix='_l', how='inner'); j = j[(j.n_e >= 10) & (j.n_l >= 10)]
    j['d'] = j.r_l / j.r_e - 1; return j
rows = []
for f in sorted(glob.glob('le/*_hist245_*.csv.gz')):
    mdl, sc, mem = os.path.basename(f)[:-7].split('_'); g = pd.read_csv(f)
    g['cell'] = CS[g.cell.values]                 # e3 writes row indices of the cells file it was given (code/cells_sel.csv)
    g = g[g.cell.isin(SEL)]
    p = periods(g, (1979, 2000), (2001, 2022))
    rows.append(dict(model=mdl, member=mem, d_med=p.d.median() * 100, dRH=(p.H_l - p.H_e).median(), dT=(p.T_l - p.T_e).median(), cells=len(p)))
R = pd.DataFrame(rows); R.to_csv('res/le_members.csv', index=False)
print(R.groupby('model').agg(n=('d_med', 'size'), d_min=('d_med', 'min'), d_med=('d_med', 'median'), d_max=('d_med', 'max'),
                             sd=('d_med', 'std'), dRH_min=('dRH', 'min'), dRH_med=('dRH', 'median'), dT=('dT', 'median')).round(2))
print('all', len(R), 'min', R.d_med.min().round(2), 'max', R.d_med.max().round(2), 'below -8.5:', (R.d_med <= -8.5).sum(), 'below -6:', (R.d_med <= -6).sum(),
      'min dRH', R.dRH.min().round(2), 'r(d,dRH)', np.corrcoef(R.d_med, R.dRH)[0, 1].round(2))
