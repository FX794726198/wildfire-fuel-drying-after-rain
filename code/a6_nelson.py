# Physics check with the Nelson (2000) dead-fuel-moisture model (hourly ERA5, NFDRS4 code; outputs of n2_nelson.py).
# For every FFMC rain event (events_obs), count days after the rain until the 1-h or 10-h stick at 14:00 local time falls to a
# dryness threshold again (same counting as the FFMC: first day k >= 1 that is dry, 31 if none within 30 days, capped at 15),
# then compare 1979-2000 with 2001-2022 per cell.
# Usage: python a6_nelson.py <dir with nelson_*.npz>
import numpy as np, pandas as pd, glob, sys
D = sys.argv[1]
EV = pd.read_csv('out_era5/events_obs.csv.gz'); SEL = pd.read_csv('cells_sel.csv').cell.values
EV = EV[EV.cell.isin(SEL)]
TH = {'m1_14': [0.08, 0.10, 0.12], 'm10_14': [0.10, 0.12, 0.15]}
rows = []
for f in sorted(glob.glob(f'{D}/nelson_*.npz')):
    z = np.load(f); cells = z['cells']; d0 = int(z['day0'])                     # column j = day index d0 + j since 1979-01-01
    for v, ths in TH.items():
        M = z[v]; nd = M.shape[1]
        for ci, c in enumerate(cells):
            e = EV[EV.cell == c]
            if e.empty: continue
            i = e.i.values; ok = (i - d0 - 1 >= 0) & (i - d0 + 31 < nd); e, i = e[ok], i[ok]
            W = M[ci][(i - d0)[:, None] + np.arange(1, 31)[None, :]]            # days i+1 .. i+30
            pre = M[ci][i - d0 - 1]                                          # day before the rain
            for t in ths:
                dry = W <= t; rec = np.where(dry.any(1), dry.argmax(1) + 1, 31)
                rows.append(pd.DataFrame(dict(cell=c, y=e.y.values, var=v, thr=t, rec=rec, pre_dry=pre <= t, rec_ffmc=e.rec.values)))
R = pd.concat(rows); R['rc'] = R.rec.clip(upper=15); R['rcf'] = R.rec_ffmc.clip(upper=15)
R.to_csv('res/nelson_events.csv.gz', index=False)
out = []
for (v, t), g in R.groupby(['var', 'thr']):
    for pre in (False, True):
        h = g[g.pre_dry] if pre else g
        k = h.assign(per=np.where(h.y < 2001, 'e', 'l')).groupby(['cell', 'per']).agg(r=('rc', 'mean'), n=('rc', 'size'), rf=('rcf', 'mean')).unstack()
        k = k[(k[('n', 'e')] >= 10) & (k[('n', 'l')] >= 10)]
        d = k[('r', 'l')] / k[('r', 'e')] - 1; df = k[('rf', 'l')] / k[('rf', 'e')] - 1
        out.append(dict(var=v, thr=t, pre_dry_only=pre, cells=len(k), d_med=d.median() * 100, shorter=(d < 0).mean(),
                        d_ffmc_same=df.median() * 100, r_cells=np.corrcoef(d.fillna(0), df.fillna(0))[0, 1],
                        mean_days_e=k[('r', 'e')].median(), r_events=np.corrcoef(h.rc, h.rcf)[0, 1]))
O = pd.DataFrame(out); O.to_csv('res/nelson_summary.csv', index=False); print(O.to_string())
