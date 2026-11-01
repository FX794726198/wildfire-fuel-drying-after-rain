# Day alignment of the MOD14A1/MYD14A1 fire pixels (classes 8-9, both sensors) with the MODIS active fire counts of the main analysis
# (dbo/fire_modis.npz from code/e4_fire.py), by correlation at shifts of -1, 0 and +1 day over all cell-days with fire in either record
# (Methods: correlation 0.91 on the same day). Run from the package root: python cloudfree/check_day_alignment.py
import numpy as np, pandas as pd, json
B = ''; pc = pd.read_csv(B + 'cloudfree/cells_cloudfree.csv').cell.values
Z = np.load(B + 'cloudfree/fm.npz'); f = Z['f89'].sum(-1)                          # [cell, day from 2002-12-01]
N = np.load(B + 'dbo/fire_modis.npz'); y0 = int(N['y0']); n = N['N'][pc]
o = (pd.Timestamp('2002-12-01') - pd.Timestamp(f'{y0}-01-01')).days; ND = f.shape[1]; n = n[:, o:o + ND]
out = {}
for L in (-1, 0, 1):                                                          # MOD14A1 day d compared with FIRMS day d + L
    a = f[:, max(0, -L):ND - max(0, L)]; b = n[:, max(0, L):ND - max(0, -L)]
    k = (a + b) > 0; out[L] = float(np.corrcoef(np.log1p(a[k]), np.log1p(b[k]))[0, 1])
best = max(out, key=out.get); print('correlation by lag', out, 'best', best)
json.dump(dict(corr=out, best=best), open(B + 'cloudfree/res/day_alignment.json', 'w'))
