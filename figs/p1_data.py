# Data for Fig. 1 (global loss of drying time after rain) and the robustness panel.
import numpy as np, pandas as pd, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import block_ci, trend_test
B = './'; O = B + 'out_era5/'; F = B + 'figs/data/'; os.makedirs(F, exist_ok=True)
C = pd.read_csv(B + 'cells_all.csv'); SEL = set(pd.read_csv(B + 'cells_sel.csv').cell)
def cells(tag, sel=True):
    R = pd.read_csv(O + f'cells_{tag}.csv').merge(C, on='cell')
    if sel: R = R[R.cell.isin(SEL)]
    R['rel'] = 100 * R.d; R['dRH'] = R.H_l - R.H_e; R['dT'] = R.T_l - R.T_e; return R
R = cells('obs'); R.to_csv(F + 'fig1_cells.csv', index=False)
print('cells', len(R), 'median', R.rel.median(), 'shorter', (R.rel < 0).mean(), 'area-weighted mean', np.average(R.rel, weights=R.area))
rows = []
for r, g in R.groupby('ar6'):
    if len(g) < 15: continue
    lo, hi = block_ci(g, nb=1000)
    rows.append(dict(ar6=r, n=len(g), med=g.rel.median(), lo=lo, hi=hi, shorter=(g.rel < 0).mean(), dRH=g.dRH.median(), dT=g.dT.median(), lat=g.lat.mean()))
G = pd.DataFrame(rows).sort_values('med'); G.to_csv(F + 'fig1_regions.csv', index=False); print(G.round(2).to_string())
# yearly global series: event drying time (capped at 15) relative to each cell's mean, averaged over cells
def series(ev, keep):
    E = ev[ev.cell.isin(keep)].copy(); E['r'] = E.rec.clip(upper=15)
    E['a'] = E.r / E.groupby('cell').r.transform('mean') - 1
    return 100 * E.groupby(['cell', 'y']).a.mean().groupby('y').mean()
S = pd.DataFrame({'era5': series(pd.read_csv(O + 'events_obs.csv.gz'), set(R.cell))})
Rc = cells('rain_chirps', sel=False); S['chirps'] = series(pd.read_csv(O + 'events_rain_chirps.csv.gz'), set(Rc.cell) & SEL)
S.to_csv(F + 'fig1_series.csv'); print({k: trend_test(S[k]) for k in S})
# robustness: medians with 5-degree block-bootstrap 95% intervals
V = [('Reference (ERA5 weather and rain)', R), ('Satellite–gauge rain (CHIRPS)', Rc[Rc.cell.isin(SEL)]), ('Gauge rain (CPC)', cells('rain_cpc')), ('Observed noon wind', cells('windobs')),
     ('Daily maximum T, minimum RH', cells('txrn')), ('Dry threshold FFMC 90', cells('rt90'))]
rows = []
for name, g in V:
    lo, hi = block_ci(g, nb=1000); rows.append(dict(version=name, n=len(g), med=g.rel.median(), lo=lo, hi=hi, shorter=(g.rel < 0).mean()))
NE = pd.read_csv(B + 'res/nelson_events.csv.gz')
for v, t, name in [('m1_14', 0.10, 'Nelson model, 1-h stick'), ('m10_14', 0.12, 'Nelson model, 10-h stick')]:
    h = NE[(NE['var'] == v) & (NE.thr == t)].assign(per=lambda x: np.where(x.y < 2001, 'e', 'l'))
    k = h.groupby(['cell', 'per']).agg(r=('rc', 'mean'), n=('rc', 'size')).unstack()
    k = k[(k[('n', 'e')] >= 10) & (k[('n', 'l')] >= 10)]
    g = pd.DataFrame(dict(cell=k.index, rel=100 * (k[('r', 'l')] / k[('r', 'e')] - 1).values)).merge(C, on='cell')
    lo, hi = block_ci(g, nb=1000); rows.append(dict(version=name, n=len(g), med=g.rel.median(), lo=lo, hi=hi, shorter=(g.rel < 0).mean()))
ST = pd.read_csv('./res/robustness_versions.csv')
for v, name in [('ERA5', 'Weather stations, ERA5 humidity'), ('Station humidity (ISD)', 'Weather stations, measured humidity')]:
    s = ST[ST.version == v].iloc[0]; rows.append(dict(version=name, n=int(s.n), med=s.med, lo=s.lo, hi=s.hi, shorter=s.shorter))
RB = pd.DataFrame(rows); RB.to_csv(F + 'fig1_robust.csv', index=False); print(RB.round(2).to_string())
