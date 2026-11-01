# Global fire response to wet fine fuels and to the wet-fuel days of each rain (ERA5 grid, MODIS/VIIRS active fires, GFED4.1s).
# A. Daily panel: Mantel-Haenszel fire-rate ratio wet (FFMC < 85) vs dry days within cell x year x month strata (fixed-effects design).
# B. Event level: fires in days 1-10 after each rain vs matched dry days of the same cell, month and year, in fire-season-mean days;
#    within-cell dose-response on wet-fuel days; implied change = slope x observed change in wet-fuel days x events per season,
#    converted to GFED carbon emissions and burned area.
# C. Gradient boosting (LightGBM): daily fire counts predicted from weather, season and place WITHOUT the wet-fuel state; the
#    observed/predicted ratio on wet vs dry days checks that A is not confounded by the covariates.
# Usage: python a4_fire.py <fire npz> <tag> [gfed npz]
import numpy as np, pandas as pd, sys, os
O = 'out_era5/'; fire_f, tag = sys.argv[1], sys.argv[2]; gfed_f = sys.argv[3] if len(sys.argv) > 3 else None
Z = np.load(fire_f); N = Z['N']; y0 = int(Z['y0']); C = pd.read_csv('cells_all.csv')
FZ = np.load(O + 'ffmc_obs.npz'); F = FZ['F'].astype(np.float32); smon = FZ['smon']
d79 = pd.date_range('1979-01-01', '2022-12-31'); off = (pd.Timestamp(f'{y0}-01-01') - d79[0]).days
nd = len(d79) - off; Nf = N[:, :nd]                                            # fire days aligned with FFMC days off..end
mon = d79.month.values[off:]; yr = d79.year.values[off:]
R = pd.read_csv(O + 'cells_obs.csv'); cells = R.cell.values
EV = pd.read_csv(O + 'events_obs.csv.gz'); EV = EV[EV.cell.isin(cells) & (EV.i - off >= 0) & (EV.i - off + 11 < nd)]
P = np.concatenate([np.load(f'era5g/era5_{y}.npz')['P24'] for y in range(1979, 2023)], 1)[:, off:]
Fw = F[:, off:]
out = {}
# ---- A. Mantel-Haenszel rate ratio within cell-year-month strata
rows = []
for c in cells:
    s = smon[c, mon - 1] & np.isfinite(Fw[c])
    if s.sum() < 100: continue
    df = pd.DataFrame(dict(y=yr[s], m=mon[s], wet=Fw[c, s] < 85, n=Nf[c, s]))
    g = df.groupby(['y', 'm', 'wet']).n.agg(['sum', 'size']).unstack('wet').dropna()
    if g.empty: continue
    a, T1, b, T0 = g[('sum', True)], g[('size', True)], g[('sum', False)], g[('size', False)]; T = T0 + T1
    rows.append(dict(cell=c, num=(a * T0 / T).sum(), den=(b * T1 / T).sum(), fires=df.n.sum()))
A = pd.DataFrame(rows).merge(C[['cell', 'ar6']], on='cell')
A = A[A.fires >= 20]
out['MH_RR_all'] = A.num.sum() / A.den.sum(); out['MH_cells'] = len(A)
out['MH_RR_cell_median'] = float(np.median((A.num / A.den).replace([np.inf], np.nan).dropna()))
A.groupby('ar6').apply(lambda x: x.num.sum() / x.den.sum()).rename('RR').to_csv(f'res/fire_MH_{tag}_regions.csv')
A.assign(RR=A.num / A.den).to_csv(f'res/fire_MH_{tag}_cells.csv', index=False)
# ---- B. event-level suppression vs matched dry days
fs_mean = {}
res = []
for c, e in EV.groupby('cell'):
    s = smon[c, mon - 1]; fire_c = Nf[c]
    # fire-season mean daily fires by year
    ym = pd.Series(fire_c[s]).groupby(yr[s]).mean()
    if ym.sum() == 0: continue
    # dry control days: < 1 mm, previous FFMC >= 85, no >= 10 mm rain in next 10 days, fire season
    cand = np.where(s & (P[c] < 1))[0]; cand = cand[(cand >= 1) & (cand + 11 < nd)]
    ok = (Fw[c, cand - 1] >= 85) & np.array([P[c, j + 1:j + 11].max() < 10 for j in cand])
    cand = cand[ok]
    W = np.array([fire_c[j + 1:j + 11].sum() for j in cand]); cy, cm = yr[cand], mon[cand]
    ctrl = pd.DataFrame(dict(y=cy, m=cm, w=W)).groupby(['y', 'm']).w.agg(['mean', 'size'])
    for _, r in e.iterrows():
        j = int(r.i) - off; key = (yr[j], mon[j])
        if key not in ctrl.index or ctrl.loc[key, 'size'] < 3: continue
        mu = ym.get(yr[j], np.nan)
        if not mu > 0: continue
        fe = fire_c[j + 1:j + 11].sum()
        res.append((c, yr[j], mon[j], r.rec, r.further10, (ctrl.loc[key, 'mean'] - fe) / mu, fe / mu, ctrl.loc[key, 'mean'] / mu))
B = pd.DataFrame(res, columns=['cell', 'y', 'm', 'rec', 'further10', 'avoid', 'fev', 'fctl']).merge(C[['cell', 'ar6']], on='cell')
B['rc'] = B.rec.clip(upper=15)
B.to_csv(f'res/fire_events_{tag}.csv.gz', index=False)
def slope(b):
    x = b.rc - b.groupby('cell').rc.transform('mean'); y = b.avoid - b.groupby('cell').avoid.transform('mean')
    return float((x * y).sum() / (x * x).sum())
out['avoid_mean'] = B.avoid.mean(); out['avoid_ratio'] = B.fev.sum() / B.fctl.sum(); out['slope'] = slope(B)
out['slope_nofurther'] = slope(B[~B.further10]); out['avoid_nofurther'] = B[~B.further10].avoid.mean()
# bootstrap over cells
rng = np.random.default_rng(1); cl = B.cell.unique(); G = {c: g for c, g in B.groupby('cell')}; bs = []
for k in range(300):
    bb = pd.concat([G[c] for c in rng.choice(cl, len(cl))]); bs.append((slope(bb), bb.avoid.mean()))
bs = np.array(bs); out['slope_ci'] = np.percentile(bs[:, 0], [2.5, 97.5]).tolist(); out['avoid_ci'] = np.percentile(bs[:, 1], [2.5, 97.5]).tolist()
# implied change per cell: slope x delta rec (capped) x events per season
Rr = R.assign(drec=R.rec15_l - R.rec15_e)
B_reg = B.groupby('ar6').apply(lambda b: pd.Series(dict(n=len(b), avoid=b.avoid.mean(), slope=slope(b) if b.cell.nunique() > 5 else np.nan)))
B_reg.to_csv(f'res/fire_dose_{tag}_regions.csv')
Rr['implied_fd'] = out['slope'] * Rr.drec * Rr.nev_l                          # fire-season-mean days of fire per season
if gfed_f:
    G5 = np.load(gfed_f); EC = G5['EC']; BA = G5['BA']; g0 = pd.Timestamp(f"{int(G5['y0'])}-01-01")
    dd = pd.date_range(g0, periods=EC.shape[1]); s0 = d79.get_loc(g0)
    mg = d79.month.values[s0:s0 + EC.shape[1]] if s0 + EC.shape[1] <= len(d79) else None
    k = min(EC.shape[1], len(d79) - s0); ECs, BAs = EC[:, :k], BA[:, :k]; mg = d79.month.values[s0:s0 + k]; yg = d79.year.values[s0:s0 + k]
    seasEC = np.array([ECs[c, smon[c, mg - 1]].mean() if smon[c].any() else 0 for c in range(len(C))])          # g C per fire-season day
    okb = yg <= 2016
    seasBA = np.array([BAs[c, smon[c, mg - 1] & okb].mean() if smon[c].any() else 0 for c in range(len(C))])    # km2 per fire-season day
    Rr['implied_TgC'] = Rr.implied_fd * seasEC[Rr.cell] / 1e12; Rr['implied_km2'] = Rr.implied_fd * seasBA[Rr.cell]
    out['implied_TgC_per_yr'] = float(Rr.implied_TgC.sum()); out['implied_Mha_per_yr'] = float(Rr.implied_km2.sum() / 1e4)
    out['seasonEC_total_TgC'] = float((seasEC[Rr.cell] * Rr.nseas.values * 30.4).sum() / 1e12)
Rr.to_csv(f'res/fire_implied_{tag}.csv', index=False)
# ---- C. LightGBM check
try:
    import lightgbm as lgb
    T12 = np.concatenate([np.load(f'era5g/era5_{y}.npz')['T12'] for y in range(1979, 2023)], 1)[:, off:]
    rs = np.random.default_rng(0); sel = rs.choice(A.cell.values, min(600, len(A)), replace=False); X = []
    for c in sel:
        s = np.where(smon[c, mon - 1] & np.isfinite(Fw[c]))[0]
        X.append(pd.DataFrame(dict(cell=c, y=yr[s], m=mon[s], T=T12[c, s], P=P[c, s], P1=P[c, np.maximum(s - 1, 0)],
                                   lat=C.lat[c], lon=C.lon[c], fy=pd.Series(Nf[c, s]).groupby(yr[s]).transform('mean').values,
                                   wet=Fw[c, s] < 85, n=Nf[c, s])))
    X = pd.concat(X); feats = ['m', 'T', 'P', 'P1', 'lat', 'lon', 'fy', 'y']
    tr = X.y % 2 == 0
    mdl = lgb.LGBMRegressor(objective='poisson', n_estimators=400, learning_rate=0.05, num_leaves=63, verbose=-1).fit(X[tr][feats], X[tr].n)
    te = X[~tr].assign(pred=mdl.predict(X[~tr][feats]))
    out['gbm_obs_over_pred_wet'] = te[te.wet].n.sum() / te[te.wet].pred.sum(); out['gbm_obs_over_pred_dry'] = te[~te.wet].n.sum() / te[~te.wet].pred.sum()
    out['gbm_ratio'] = out['gbm_obs_over_pred_wet'] / out['gbm_obs_over_pred_dry']
except Exception as ex: out['gbm_error'] = str(ex)
pd.Series(out).to_json(f'res/fire_summary_{tag}.json', indent=1); print(tag, out)
