# Fig. 3b and Fig. S4b-d on the fixed-event design (v1.2): changes after the same rain events with trends removed (fx/ runs).
import numpy as np, pandas as pd, os
B = './'; F = B + 'figs/data/'; R = B + 'res/'
C = pd.read_csv(B + 'code/cells_sel.csv').set_index('cell')
def percell(ref, v, cap=31):
    ev = pd.read_csv(f'{B}fx/events_{ref}.csv.gz'); r = np.minimum(np.load(f'{B}fx/rec_{ref}_{v}.npy').astype(float), cap)
    g = pd.DataFrame(dict(cell=ev.cell.values, late=ev.y.values >= 2001, r=r)).groupby(['cell', 'late']).r.agg(['mean', 'size']).unstack()
    g = g[(g[('size', False)] >= 10) & (g[('size', True)] >= 10) & (g[('size', False)] + g[('size', True)] >= 30)]
    return (g[('mean', True)] / g[('mean', False)] - 1) * 100
rows = []; lab = {'obs': 'Observed', 'warm': 'Warming only', 'dry': 'Drying only', 'neither': 'Neither'}
for rain, ref in [('ERA5 rain', 'noon'), ('CHIRPS rain', 'chirps')]:
    S = pd.read_csv(f'{R}fixed_{ref}_cap31.csv').set_index('variant')
    for v in ['obs', 'warm', 'dry', 'neither']:
        s = S.loc[v]; rows.append(dict(rain=rain, case=lab[v], n=int(s.cells), med=s.med, lo=s.lo, hi=s.hi, paired=s.paired, p_lo=s.p_lo, p_hi=s.p_hi, shorter=s.shorter))
DEC = pd.DataFrame(rows); DEC.to_csv(F + 'fig3_decomp.csv', index=False); print(DEC.round(2).to_string())
M = pd.read_csv(F + 'figS3_maps.csv').set_index('cell') if 'cell' in pd.read_csv(F + 'figS3_maps.csv', nrows=1) else None
X = pd.DataFrame({v: percell('noon', v) for v in ['warm', 'dry', 'neither']}).join(C[['lat', 'lon']])
old = pd.read_csv(F + 'figS3_maps.csv')
X = X.reset_index().rename(columns={'index': 'cell'}).merge(old[['cell', 'dT']], on='cell', how='left') if 'cell' in old else X
X.to_csv(F + 'figS3_maps.csv', index=False); print(X.describe().round(2))
