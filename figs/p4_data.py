# Data for Fig. 4 and Figs. S6-S7 (human influence, climate models, projections) from the CMIP6 summary of a5_cmip6.py.
# Usage: python p4_data.py [tag]   (tag of the res/<tag>_*.csv files; 'test' = CanESM5 only, 'cmip6' = all models)
import numpy as np, pandas as pd, sys, os, json
B = './'; R = B + 'res/'; F = B + 'figs/data/'; tag = sys.argv[1] if len(sys.argv) > 1 else 'cmip6'
C = pd.read_csv(B + 'cells_all.csv')
H = pd.read_csv(R + f'{tag}_hist.csv')
OBS = pd.read_csv(R + 'test_hist.csv').query("model == 'ERA5'").iloc[0]            # ERA5 daily max T + min RH, same counting as the models
Hm = H[H.model != 'ERA5'].copy()
# one row per model for the comparison with observations: SSP2-4.5 extension where available, otherwise SSP5-8.5 (they differ only after 2014)
Hm['use'] = False; Hm.loc[Hm.sort_values('scen').groupby('model').head(1).index, 'use'] = True
Hm.to_csv(F + 'fig4_hist.csv', index=False); Hu = Hm[Hm.use]
P = pd.read_csv(R + f'{tag}_proj.csv')
gcols = [c for c in P.columns if c.startswith('gwl_')]
# one value per model and warming level: mean over the scenarios that reach it
G = P.groupby('model')[gcols].mean()
G.to_csv(F + 'fig4_gwl_models.csv')
pi = R + f'{tag}_pi.csv'
if os.path.exists(pi):
    Pi = pd.read_csv(pi); Pi.to_csv(F + 'fig4_pi.csv', index=False)
else:
    Pi = None
    if os.path.exists(F + 'fig4_pi.csv'): os.remove(F + 'fig4_pi.csv')
# per-cell multi-model median at each warming level (map) and model agreement
PC = pd.read_csv(R + f'{tag}_proj_cells.csv.gz')
PC['d'] = 100 * PC.d
M = PC.groupby(['cell', 'model', 'win']).d.mean().reset_index()                     # average the two scenarios
MM = M.groupby(['cell', 'win']).agg(med=('d', 'median'), n=('d', 'size'), neg=('d', lambda x: (x < 0).mean())).reset_index()
MM = MM.merge(C[['cell', 'lat', 'lon', 'ar6']], on='cell'); MM.to_csv(F + 'fig4_map.csv', index=False)
S = pd.read_csv(R + f'{tag}_series.csv', index_col=0)
# observed series with the same counting as the models (ERA5 daily max T + min RH; per-cell relative anomaly of event-weighted rec15)
SEL = set(pd.read_csv(B + 'cells_sel.csv').cell)
E = pd.read_csv(B + 'out_era5/events_txrn.csv.gz'); E = E[E.cell.isin(SEL)]
og = E.groupby(['cell', 'y']).agg(rec15=('rec', lambda x: x.clip(upper=15).mean()), n=('rec', 'size')).reset_index()
mref = og[(og.y >= 1979) & (og.y <= 2022)].groupby('cell').rec15.mean(); og['a'] = og.rec15 / og.cell.map(mref) - 1
S['ERA5_obs'] = og.groupby('y').a.mean(); S.to_csv(F + 'figS6_series.csv')
CE = pd.read_csv(R + f'{tag}_cells.csv.gz'); CE['d'] = 100 * CE.d; CE.to_csv(F + 'figS6_cells.csv.gz', index=False)
summ = dict(obs_d=OBS.d_med, obs_dRH=OBS.dRH, obs_dT=OBS.dT, obs_trend=OBS.trend, n_models=int(Hm.model.nunique()),
            hist_d_med=Hu.d_med.median(), hist_d_min=Hu.d_med.min(), hist_d_max=Hu.d_med.max(),
            hist_dRH_med=Hu.dRH.median(), hist_dRH_min=Hu.dRH.min(), hist_dRH_max=Hu.dRH.max(), hist_dT_med=Hu.dT.median(), hist_dT_min=Hu.dT.min(), hist_dT_max=Hu.dT.max(),
            hist_trend_med=Hu.trend.median(), r_d_dRH=float(np.corrcoef(Hu.dRH, Hu.d_med)[0, 1]), r_d_dT=float(np.corrcoef(Hu.dT, Hu.d_med)[0, 1]),
            gwl={c: [G[c].median(), G[c].quantile(0.1), G[c].quantile(0.9), int(G[c].notna().sum())] for c in gcols})
if Pi is not None:
    v = Pi.d_med.values; summ.update(pi_sd=float(np.std(v)), pi_min=float(v.min()), pi_p_obs=float((v <= OBS.d_med).mean()), pi_p_models=float((v <= summ['hist_d_med']).mean()),
                                    pi_sd_trend=float(Pi.trend.std()), pi_trend_min=float(Pi.trend.min()),
                                    pi_q=[float(q) for q in np.percentile(v, [1, 2.5, 5, 95, 97.5, 99])], pi_models=int(Pi.model.nunique()), pi_windows=len(v))
json.dump(summ, open(F + 'fig4_summary.json', 'w'), indent=1, default=float); print(json.dumps(summ, indent=1, default=float))
