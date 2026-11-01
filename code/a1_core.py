# Core global numbers from the ERA5 grid run: time series + trend, rain response and compensation factor, analytic scaling.
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
O = 'out_era5/'; C = pd.read_csv('cells_all.csv')
R = pd.read_csv(O + 'cells_obs.csv').merge(C, on='cell'); E = pd.read_csv(O + 'events_obs.csv.gz'); E = E[E.cell.isin(R.cell)]
def ar1_trend(y, x):
    b, a, r, p, se = stats.linregress(x, y); res = y - (a + b * x); r1 = np.corrcoef(res[:-1], res[1:])[0, 1]
    neff = len(y) * (1 - r1) / (1 + r1); t = b / se * np.sqrt(max(neff - 2, 1) / (len(y) - 2)); return b, 2 * stats.t.sf(abs(t), max(neff - 2, 1))
E['a'] = E.rec.clip(upper=15) / E.groupby('cell').rec.transform(lambda x: x.clip(upper=15).mean()) - 1
s = E.groupby(['cell', 'y']).a.mean().groupby('y').mean()
b, p = ar1_trend(s.values, s.index.values.astype(float)); print('global yearly rel anomaly trend %/dec', round(1000 * b, 2), 'p', p)
s.to_csv('res/series_global.csv')
# rain response: rec ~ ln p + late + cell FE
E['lp'] = np.log(E.p); E['late'] = (E.y >= 2001).astype(float); E['r15'] = E.rec.clip(upper=15)
X = E[['lp', 'late']] - E.groupby('cell')[['lp', 'late']].transform('mean'); Y = E.r15 - E.groupby('cell').r15.transform('mean')
m = sm.OLS(Y, X).fit(cov_type='cluster', cov_kwds={'groups': E.cell}); print(m.params, m.bse)
beta, gam = m.params.lp, m.params.late; print('rain factor to compensate', np.exp(-gam / beta))
print('CC 7%/K x 0.78 K gives rain factor', round(1.07 ** 0.78, 3), '-> days', round(beta * np.log(1.07 ** 0.78), 3), 'vs loss', round(gam, 3))
# weights
for w in ['area', 'pop2020']:
    print(w, 'weighted mean d', round(np.average(R.d, weights=R[w] + 1e-9), 4))
print('share of cells shorter', (R.d < 0).mean(), 'share of area', R.area[R.d < 0].sum() / R.area.sum(), 'pop living in cells shorter (bn)', R.pop2020[R.d < 0].sum() / 1e9, 'of', R.pop2020.sum() / 1e9)
