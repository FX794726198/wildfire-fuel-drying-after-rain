# Per-model diagnostics for the CMIP6 comparison (Table S5): with each model's own fire season (v1.1 runs, res/cmip6_hist.csv) and
# with the observed fire-season months (fs/ runs): cells with >= 10 events in each period, events, baseline (1979-2000) mean drying
# time (capped at 15 d), change in all of the model's qualifying cells and in the cells common to ERA5 and all models, and the
# fire-season changes in daily min RH and max T in the common cells.  ERA5 (daily max T / min RH) in the first row.
import numpy as np, pandas as pd, glob, os, sys
D = '.'; CS = pd.read_csv(f'{D}/code/cells_sel.csv'); SEL = CS.cell.values
common = pd.read_csv(f'{D}/res/samecov_common.csv', index_col=0).index
def stats(cell, y, s, n, cells=None):
    t = pd.DataFrame(dict(cell=cell, late=y >= 2001, s=s, n=n)).groupby(['cell', 'late'])[['s', 'n']].sum().unstack()
    t = t[(t[('n', False)] >= 10) & (t[('n', True)] >= 10)]
    d = (t[('s', True)] / t[('n', True)]) / (t[('s', False)] / t[('n', False)]) * 100 - 100
    base = (t[('s', False)] / t[('n', False)])
    return dict(cells=len(t), events=int(t[('n', False)].sum() + t[('n', True)].sum()), base=base.median(), d_all=d.median(),
                d_common=d.reindex(common).median(), base_common=base.reindex(common).median())   # baseline in the common cells (v1.3)
rows = []
ev = pd.read_csv(f'{D}/fx/events_txrn.csv.gz'); r = np.minimum(np.load(f'{D}/fx/rec_txrn_obs.npy').astype(float), 15)
rr = dict(model='ERA5', **stats(ev.cell.values, ev.y.values, r, np.ones(len(r))))
# ERA5 fire-season daily min RH / max T change in the common cells (observed season months)
smon = np.load(f'{D}/out_era5/ffmc_obs.npz')['smon'][SEL]; dates = pd.date_range('1979-01-01', '2022-12-31'); mon, yrs = dates.month.values, dates.year.values
Z = [np.load(f'{D}/era5g/era5_{y}.npz') for y in range(1979, 2023)]
row = {c: j for j, c in enumerate(SEL)}; ic = np.array([row[c] for c in common])
Tx = np.concatenate([z['Tmax'][SEL[ic]] for z in Z], 1); Hn = np.concatenate([z['RHmin'][SEL[ic]] for z in Z], 1); del Z
seas = smon[ic][:, mon - 1]
chg = lambda X: np.nanmean(np.where(seas & (yrs >= 2001), X, np.nan), 1) - np.nanmean(np.where(seas & (yrs < 2001), X, np.nan), 1)
rr.update(dRH=np.median(chg(Hn)), dT=np.median(chg(Tx))); rows.append(rr)
OWN = pd.read_csv(f'{D}/res/cmip6_hist.csv')
for f in sorted(glob.glob(f'{D}/fs/*_hist*_*.csv.gz')):
    mdl, scen, mem = os.path.basename(f)[:-7].split('_'); g = pd.read_csv(f); g = g[(g.y >= 1979) & (g.y <= 2022)]
    ge = g[g.n.notna() & (g.n > 0)]
    st = stats(SEL[ge.cell.values], ge.y.values, (ge.rec15 * ge.n).values, ge.n.values)
    gc = g[np.isin(SEL[g.cell.values], common)]
    hc = gc.groupby(['cell', gc.y >= 2001])[['H', 'T']].mean().unstack()
    o = OWN[(OWN.model == mdl) & (OWN.scen == scen)]
    rows.append(dict(model=mdl, scen=scen, member=mem, **st, dRH=(hc[('H', True)] - hc[('H', False)]).median(), dT=(hc[('T', True)] - hc[('T', False)]).median(),
                     d_own=o.d_med.iloc[0] if len(o) else np.nan, cells_own=o.cells.iloc[0] if len(o) else np.nan))
T = pd.DataFrame(rows); T.to_csv(f'{D}/res/models_table.csv', index=False); print(T.round(2).to_string())
