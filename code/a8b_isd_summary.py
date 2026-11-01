# Summarise the station-vs-ERA5 comparison: medians with bootstrap intervals (resampling stations), by zone, US vs rest,
# with stricter coverage and after removing stations with a break in the station-minus-ERA5 humidity difference.
import numpy as np, pandas as pd, sys, json, os
tag = sys.argv[1] if len(sys.argv) > 1 else 'main'; MINCOV = float(os.environ.get('MINCOV', 0.7))
R = pd.read_csv(f'./res/isd_{tag}.csv'); rng = np.random.default_rng(1)
def snht(x):
    x = x[np.isfinite(x)]; n = len(x)
    if n < 20: return np.nan, np.nan
    z = (x - x.mean()) / x.std(); best, sh = 0, 0
    for k in range(5, n - 5):
        t = k * z[:k].mean() ** 2 + (n - k) * z[k:].mean() ** 2
        if t > best: best, sh = t, x[k:].mean() - x[:k].mean()
    return best, sh
dif = R[[f'dif{y}' for y in range(1979, 2023)]].values
R['snht'], R['shift'] = zip(*[snht(r) for r in dif])
R['brk'] = (R.snht > 8.6) & (R['shift'].abs() > 2)              # 95% SNHT critical value for n~44; shift > 2 pp
R['zone'] = pd.cut(R.lat, [-90, -23.5, 23.5, 90], labels=['S extratropics', 'tropics', 'N extratropics']).astype(str)
base = (R.cov_e >= MINCOV) & (R.cov_l >= MINCOV) & R.d_st.notna() & R.d_era.notna()
def ci(x, f=np.median, B=2000):
    x = np.asarray(x); b = [f(rng.choice(x, len(x))) for _ in range(B)]; return [float(f(x)), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]
out = {}
sets = {'all': base, 'cov90': base & (R.cov_e >= 0.9) & (R.cov_l >= 0.9), 'no_break': base & ~R.brk, 'non_US': base & (R.ctry != 'US'),
        'US': base & (R.ctry == 'US'), 'tropics': base & (R.zone == 'tropics'), 'N extratropics': base & (R.zone == 'N extratropics'),
        'S extratropics': base & (R.zone == 'S extratropics'), 'tropics_no_break': base & (R.zone == 'tropics') & ~R.brk}
for k, m in sets.items():
    x = R[m]
    out[k] = dict(n=int(m.sum()), d_st=ci(x.d_st * 100), d_era=ci(x.d_era * 100), diff=ci((x.d_st - x.d_era) * 100),
                  dRH_st=ci(x.dRH_st), dRH_era=ci(x.dRH_era), dRH_diff=ci(x.dRH_st - x.dRH_era), dT_st=ci(x.dT_st), dT_era=ci(x.dT_era),
                  shorter_st=float((x.d_st < 0).mean()), shorter_era=float((x.d_era < 0).mean()),
                  r_d=float(np.corrcoef(x.d_st, x.d_era)[0, 1]), r_dRH=float(np.corrcoef(x.dRH_st, x.dRH_era)[0, 1]),
                  n_brk=int(x.brk.sum()))
    print(f"{k:16s} n={out[k]['n']:4d} d_st {out[k]['d_st'][0]:6.2f} [{out[k]['d_st'][1]:6.2f},{out[k]['d_st'][2]:6.2f}]  d_era {out[k]['d_era'][0]:6.2f} [{out[k]['d_era'][1]:6.2f},{out[k]['d_era'][2]:6.2f}]  "
          f"dRH st {out[k]['dRH_st'][0]:5.2f} era {out[k]['dRH_era'][0]:5.2f} diff {out[k]['dRH_diff'][0]:5.2f} [{out[k]['dRH_diff'][1]:5.2f},{out[k]['dRH_diff'][2]:5.2f}]  dT st {out[k]['dT_st'][0]:4.2f} era {out[k]['dT_era'][0]:4.2f}  r_d {out[k]['r_d']:.2f} r_dRH {out[k]['r_dRH']:.2f} brk {out[k]['n_brk']}")
R.drop(columns=[f'dif{y}' for y in range(1979, 2023)]).to_csv(f'./res/isd_{tag}_flags.csv', index=False)
json.dump(out, open(f'./res/isd_{tag}_summary.json', 'w'), indent=1)
