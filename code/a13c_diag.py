# Across-model diagnostics for the fixed-season comparison: correlation of the change in the common cells with the simulated humidity
# and temperature changes and with the baseline drying time; regression on dRH and baseline.
import pandas as pd, numpy as np, statsmodels.api as sm
T = pd.read_csv('./res/models_table.csv'); T['base'] = T.base_common   # v1.3: baseline over the same common cells as the change
E = T[T.model == 'ERA5'].iloc[0]; M = T[T.model != 'ERA5']
print('n models', len(M), 'cells', M.cells.min(), M.cells.max(), 'events', M.events.min(), M.events.max(), 'base', M.base.min().round(2), M.base.max().round(2), 'ERA5 base', round(E.base, 2))
print('median d_all', M.d_all.median().round(2), M.d_all.min().round(2), M.d_all.max().round(2), '| d_common', M.d_common.median().round(2), M.d_common.min().round(2), M.d_common.max().round(2))
print('dRH', M.dRH.min().round(2), M.dRH.max().round(2), 'dT', M.dT.min().round(2), M.dT.max().round(2), 'ERA5', round(E.dRH, 2), round(E.dT, 2))
for c in ['dRH', 'dT', 'base']: print('r(d_common,', c, ')', round(np.corrcoef(M.d_common, M[c])[0, 1], 2))
X = sm.add_constant(M[['dRH', 'base']]); f = sm.OLS(M.d_common, X).fit(); print(f.params.round(2).to_dict(), f.pvalues.round(3).to_dict())
pred = f.params['const'] + f.params['dRH'] * E.dRH + f.params['base'] * E.base; print('predicted at ERA5 dRH and base', round(pred, 2))
s1 = np.polyfit(M.dRH, M.d_common, 1); print('dRH-only fit at ERA5 dRH', round(np.polyval(s1, E.dRH), 2))
print('own-season d_own median', M.d_own.median().round(2), 'fixed vs own r', round(np.corrcoef(M.d_own, M.d_all)[0, 1], 2))
import json
json.dump(dict(n=len(M), r_dRH=float(np.corrcoef(M.d_common, M.dRH)[0, 1]), r_dT=float(np.corrcoef(M.d_common, M.dT)[0, 1]),
               r_base=float(np.corrcoef(M.d_common, M.base)[0, 1]), pred=float(pred), dRH_min=float(M.dRH.min()), dRH_max=float(M.dRH.max()),
               dT_min=float(M.dT.min()), dT_max=float(M.dT.max()), era_dRH=float(E.dRH), era_dT=float(E.dT), base_min=float(M.base.min()),
               base_max=float(M.base.max()), era_base=float(E.base), own_med=float(M.d_own.median()), base_coef=float(f.params['base']),
               base_p=float(f.pvalues['base']), base_t=float(f.tvalues['base']), r_own=float(np.corrcoef(M.d_own, M.d_all)[0, 1])), open('./res/models_diag.json', 'w'), indent=1)
