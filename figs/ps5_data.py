# Data for the station humidity check (new Fig. S5, Fig. 1b rows, Table S1): station vs ERA5 at the same cell and hour.
# Medians with 5-degree block-bootstrap 95% intervals, as for the gridded results.
import numpy as np, pandas as pd, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import block_ci
C = './'; D = os.path.dirname(os.path.abspath(__file__)) + '/data/'
out = {}
for tag in ['main', 'cpc']:
    R = pd.read_csv(C + f'res/isd_{tag}_flags.csv')
    k = (R.cov_e >= 0.7) & (R.cov_l >= 0.7) & R.d_st.notna() & R.d_era.notna(); X = R[k].copy()
    X['st'] = X.d_st * 100; X['era'] = X.d_era * 100; X['diff'] = X.st - X.era; X['dh'] = X.dRH_st - X.dRH_era
    X.to_csv(D + f'figS5_stations_{tag}.csv', index=False)
    grp = {'all': X, 'no_break': X[~X.brk], 'non_US': X[X.ctry != 'US'], 'N extratropics': X[X.zone == 'N extratropics'],
           'S extratropics': X[X.zone == 'S extratropics'], 'tropics': X[X.zone == 'tropics'], 'cov90': X[(X.cov_e >= 0.9) & (X.cov_l >= 0.9)]}
    o = {}
    for g, Y in grp.items():
        o[g] = dict(n=len(Y), countries=int(Y.ctry.nunique()))
        for col in ['st', 'era', 'diff', 'dRH_st', 'dRH_era', 'dh', 'dT_st', 'dT_era']:
            lo, hi = block_ci(Y, col=col); o[g][col] = [float(Y[col].median()), float(lo), float(hi)]
        o[g]['shorter_st'] = float((Y.st < 0).mean()); o[g]['shorter_era'] = float((Y.era < 0).mean())
        o[g]['r_dRH'] = float(np.corrcoef(Y.dRH_st, Y.dRH_era)[0, 1]); o[g]['r_d'] = float(np.corrcoef(Y.st, Y.era)[0, 1])
    out[tag] = o
    print(tag, {g: (o[g]['n'], [round(v, 2) for v in o[g]['st']], [round(v, 2) for v in o[g]['era']], [round(v, 2) for v in o[g]['diff']], [round(v, 2) for v in o[g]['dh']]) for g in o})
json.dump(out, open(D + 'figS5_summary.json', 'w'), indent=1)
