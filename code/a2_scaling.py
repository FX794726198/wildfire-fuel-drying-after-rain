# Analytic scaling of wet-fuel days with noon RH and T from the FFMC drying equation (Van Wagner 1987):
# N ~ log10((m0 - E)/(m* - E)) / k_d, E = equilibrium moisture (drying), k_d = log drying rate; m* = moisture at FFMC 85.
import numpy as np, pandas as pd
def N(T, H, W=12., m0=60., F=85.):
    ms = 147.2 * (101 - F) / (59.5 + F)
    ed = 0.942 * H ** 0.679 + 11 * np.exp((H - 100) / 10) + 0.18 * (21.1 - T) * (1 - np.exp(-0.115 * H))
    kd = (0.424 * (1 - (H / 100) ** 1.7) + 0.0694 * np.sqrt(W) * (1 - (H / 100) ** 8)) * 0.581 * np.exp(0.0365 * T)
    return np.log10((m0 - ed) / np.maximum(ms - ed, 0.5)) / kd
C = pd.read_csv('cells_all.csv'); R = pd.read_csv('out_era5/cells_obs.csv').merge(C, on='cell')
R['pred'] = N(R.T_l, R.H_l) / N(R.T_e, R.H_e) - 1
R['predH'] = N(R.T_e, R.H_l) / N(R.T_e, R.H_e) - 1; R['predT'] = N(R.T_l, R.H_e) / N(R.T_e, R.H_e) - 1
R['sH'] = (N(R.T_e, R.H_e - 1) / N(R.T_e, R.H_e) - 1); R['sT'] = (N(R.T_e + 1, R.H_e) / N(R.T_e, R.H_e) - 1)
print('median sens %/pp drying', round(100 * R.sH.median(), 2), '%/K', round(100 * R.sT.median(), 2))
print('median predicted', round(R.pred.median(), 4), 'H part', round(R.predH.median(), 4), 'T part', round(R.predT.median(), 4), 'observed', round(R.d.median(), 4))
g = R.groupby('ar6').filter(lambda x: len(x) >= 15).groupby('ar6').agg(n=('d', 'size'), obs=('d', 'median'), pred=('pred', 'median'))
print(g.round(3).to_string()); print('r regions', np.corrcoef(g.obs, g.pred)[0, 1], 'r cells', np.corrcoef(R.d, R.pred)[0, 1])
# quantile bins of predicted change
R['q'] = pd.qcut(R.pred, 10); print(R.groupby('q').agg(obs=('d', 'median'), pred=('pred', 'median'), n=('d', 'size')).round(3))
R[['cell', 'ar6', 'd', 'pred', 'predH', 'predT', 'sH', 'sT']].to_csv('res/scaling_cells.csv', index=False); g.to_csv('res/scaling_regions.csv')
