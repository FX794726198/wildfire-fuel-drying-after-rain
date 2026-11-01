# Data for the daily max T / min RH station rows of Fig. S5d (v1.2): versions A (ERA5-assisted) and B (station-only), res/isd_txrn.csv
import numpy as np, pandas as pd, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import block_ci
C = './'; D = os.path.dirname(os.path.abspath(__file__)) + '/data/'
R = pd.read_csv(C + 'res/isd_txrn.csv'); R['zone'] = np.where(R.lat > 23.5, 'N extratropics', np.where(R.lat < -23.5, 'S extratropics', 'tropics'))
okA = R.dA_st.notna() & R.dA_era.notna() & (R.cov_e >= 0.7) & (R.cov_l >= 0.7); okB = R.dB_st.notna() & R.dB_era.notna()
out = {}
for v, ok in [('A', okA), ('B', okB)]:
    X = R[ok].copy(); X['st'] = X[f'd{v}_st'] * 100; X['era'] = X[f'd{v}_era'] * 100; X['diff'] = X.st - X.era
    if v == 'B':                                       # version B: T/RH changes from the station-only values (v1.3)
        for c in ['dRH_st', 'dRH_era', 'dT_st', 'dT_era']: X[c] = X[c.replace('_st', '_stB').replace('_era', '_eraB')]
    grp = {'all': X, 'non_US': X[X.ctry != 'US'], 'N extratropics': X[X.zone == 'N extratropics'], 'S extratropics': X[X.zone == 'S extratropics'], 'tropics': X[X.zone == 'tropics']}
    for g, Y in grp.items():
        o = dict(n=len(Y), countries=int(Y.ctry.nunique()), events=int((Y[f'ne{v}'] + Y[f'nl{v}']).sum()))
        for col in ['st', 'era', 'diff', 'dRH_st', 'dRH_era', 'dT_st', 'dT_era']:
            lo, hi = block_ci(Y, col=col); o[col] = [float(Y[col].median()), float(lo), float(hi)]
        out[f'{v}_{g}'] = o
json.dump(out, open(D + 'figS5_txrn.json', 'w'), indent=1)
for k, o in out.items(): print(k, o['n'], o['events'], [round(x, 1) for x in o['st']], [round(x, 1) for x in o['era']], [round(x, 1) for x in o['diff']], round(o['dRH_st'][0], 2), round(o['dRH_era'][0], 2))
