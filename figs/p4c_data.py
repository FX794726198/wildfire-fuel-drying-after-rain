# Data for Fig. 4a (DAMIP on fixed events, corrected single-forcing definitions) and Fig. 4d (same-coverage comparison), v1.2.
import numpy as np, pandas as pd, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from figcommon import block_ci
B = './'; R = B + 'res/'; F = B + 'figs/data/'
T = pd.read_csv(R + 'fixed_noon_cap31_dmp.csv')
T['model'] = T.variant.str.split('_').str[1]; T['kind'] = T.variant.str.split('_', n=2).str[2]
T.to_csv(F + 'fig4_damip.csv', index=False)
print(T.groupby('kind')[['med', 'paired']].agg(['median', 'min', 'max']).round(2))
S = json.load(open(R + 'samecov.json')); M = pd.read_csv(R + 'samecov_models.csv', index_col=0)
json.dump(S, open(F + 'fig4_samecov.json', 'w'), indent=1); M.to_csv(F + 'fig4_samecov_models.csv')
# v1.3: 95% interval of the ERA5 change with each cell's own fire season (the star in Fig. 4b,c,e), from resampling 5-degree blocks
E = pd.read_csv('./out_era5/events_txrn.csv.gz'); E = E[E.cell.isin(set(pd.read_csv('./code/cells_sel.csv').cell))]
E['r'] = np.minimum(E.rec, 15); E['late'] = E.y >= 2001; t = E.groupby(['cell', 'late']).r.agg(['sum', 'count']).unstack()
t = t[(t[('count', False)] >= 10) & (t[('count', True)] >= 10)]
X = pd.DataFrame({'d': (t[('sum', True)] / t[('count', True)]) / (t[('sum', False)] / t[('count', False)]) * 100 - 100}).join(pd.read_csv('./cells_all.csv').set_index('cell')[['lat', 'lon']])
json.dump(dict(obs_d=float(X.d.median()), obs_d_ci=[float(x) for x in block_ci(X, 'd', nb=1000)]), open(F + 'fig4_obs_ci.json', 'w'), indent=1)
