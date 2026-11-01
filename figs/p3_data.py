# Data for Fig. 3 (why: warmer, drier air; analytic drying law; rain cannot compensate).
import numpy as np, pandas as pd, sys, os, statsmodels.api as sm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import block_ci
B = './'; O = B + 'out_era5/'; F = B + 'figs/data/'
C = pd.read_csv(B + 'cells_all.csv'); SEL = set(pd.read_csv(B + 'cells_sel.csv').cell)
def cells(tag):
    R = pd.read_csv(O + f'cells_{tag}.csv').merge(C, on='cell'); R = R[R.cell.isin(SEL)]; R['rel'] = 100 * R.d; return R
rows = []
for rain, tags in [('ERA5 rain', ['obs', 'warm', 'dry', 'neither']), ('CHIRPS rain', ['rain_chirps', 'rain_chirps_warm', 'rain_chirps_dry', 'rain_chirps_neither'])]:
    for lab, t in zip(['Observed', 'Warming only', 'Drying only', 'Neither'], tags):
        R = cells(t); lo, hi = block_ci(R, nb=1000)
        rows.append(dict(rain=rain, case=lab, n=len(R), med=R.rel.median(), lo=lo, hi=hi, shorter=(R.rel < 0).mean(), dRH=(R.H_l - R.H_e).median(), dT=(R.T_l - R.T_e).median()))
DEC = pd.DataFrame(rows); DEC.to_csv(F + 'fig3_decomp.csv', index=False); print(DEC.round(2).to_string())
# analytic scaling deciles
SC = pd.read_csv(B + 'res/scaling_cells.csv'); SC['q'] = pd.qcut(SC.pred, 10, labels=False)
Q = SC.groupby('q').agg(pred=('pred', 'median'), obs=('d', 'median'), n=('d', 'size'))
bs = []
rng = np.random.default_rng(0)
for _ in range(1000):
    k = SC.sample(len(SC), replace=True, random_state=int(rng.integers(1e9))); bs.append(k.groupby('q').d.median().reindex(Q.index).values)
Q['lo'], Q['hi'] = np.nanpercentile(bs, 2.5, 0), np.nanpercentile(bs, 97.5, 0)
Q.to_csv(F + 'fig3_scaling.csv'); print(Q.round(3)); print('sens', SC.sH.median(), SC.sT.median(), 'r', np.corrcoef(SC.d, SC.pred)[0, 1])
# observed fire-season noon RH change per cell for the map
R = cells('obs'); R['dRH'] = R.H_l - R.H_e; R['dT'] = R.T_l - R.T_e; R[['cell', 'lat', 'lon', 'dRH', 'dT', 'rel']].to_csv(F + 'fig3_cells.csv', index=False)
# rain compensation: drying time vs rain amount, early and late (cell fixed effects removed)
E = pd.read_csv(O + 'events_obs.csv.gz'); E = E[E.cell.isin(set(R.cell))].copy()
E['r15'] = E.rec.clip(upper=15); E['late'] = (E.y >= 2001).astype(float); E['lp'] = np.log(E.p)
X = E[['lp', 'late']] - E.groupby('cell')[['lp', 'late']].transform('mean'); Y = E.r15 - E.groupby('cell').r15.transform('mean')
mo = sm.OLS(Y, X).fit(cov_type='cluster', cov_kwds={'groups': E.cell})
beta, gam = mo.params.lp, mo.params.late; fac = np.exp(-gam / beta)
# delta-method-free CI for the factor by cluster bootstrap of cells
cl = E.cell.unique(); Gc = {c: i for c, i in E.groupby('cell').groups.items()}; fs = []
for _ in range(300):
    idx = np.concatenate([Gc[c] for c in rng.choice(cl, len(cl))]); x, y = X.loc[idx], Y.loc[idx]
    b = np.linalg.lstsq(x.values, y.values, rcond=None)[0]; fs.append(np.exp(-b[1] / b[0]))
print('beta', beta, 'gamma', gam, 'factor', fac, np.percentile(fs, [2.5, 97.5]))
E['res'] = E.r15 - E.groupby('cell').r15.transform('mean') + E.r15.mean()
bins = [10, 12.5, 15, 20, 25, 35, 50, 80, 1e4]; E['pb'] = pd.cut(E.p, bins)
CU = E.groupby(['pb', 'late'], observed=True).agg(p=('p', 'median'), r=('res', 'mean'), se=('res', 'sem'), n=('res', 'size')).reset_index()
CU.to_csv(F + 'fig3_rain_curves.csv', index=False); print(CU.round(3).to_string())
pd.Series(dict(beta=beta, gamma=gam, factor=fac, f_lo=np.percentile(fs, 2.5), f_hi=np.percentile(fs, 97.5), mean_r15=E.r15.mean(),
               cc_days=beta * np.log(1.07 ** R.eval('T_l - T_e').median()), dT=R.eval('T_l - T_e').median())).to_csv(F + 'fig3_rain_fit.csv')
