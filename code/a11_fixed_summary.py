# Summaries of the fixed-event variants (e7_fixed.py): per-cell change in mean drying time between 1979-2000 and 2001-2022 after the
# same rain events, median over the same cells, and the paired difference from the reference (median over cells of the per-cell
# difference), with 95% intervals from resampling 5-degree blocks of cells.  Usage: python a11_fixed_summary.py <ref> [cap]
import numpy as np, pandas as pd, glob, os, sys, json
sys.path.insert(0, './figs'); from figcommon import block_ci
D = '.'; ref = sys.argv[1]; cap = int(sys.argv[2]) if len(sys.argv) > 2 else 31
ev = pd.read_csv(f'{D}/fx/events_{ref}.csv.gz'); C = pd.read_csv(f'{D}/code/cells_sel.csv').set_index('cell')
YEND = int(os.environ.get('YEND', 2022)); keep = ev.y.values <= YEND; ev = ev[keep]       # YEND=2020: no extension of DAMIP runs
PAT = os.environ.get('PAT', '*')
late = ev.y.values >= 2001
def percell(r):
    r = np.minimum(r[keep].astype(float), cap)
    g = pd.DataFrame(dict(cell=ev.cell.values, late=late, r=r)).groupby(['cell', 'late']).r.agg(['mean', 'size']).unstack()
    g = g[(g[('size', False)] >= 10) & (g[('size', True)] >= 10) & (g[('size', False)] + g[('size', True)] >= 30)]
    return (g[('mean', True)] / g[('mean', False)] - 1) * 100
R0 = percell(np.load(f'{D}/fx/rec_{ref}_obs.npy'))
rows = []
for f in sorted(glob.glob(f'{D}/fx/rec_{ref}_{PAT}.npy')):
    v = os.path.basename(f)[len(f'rec_{ref}_'):-4]; d = percell(np.load(f))
    j = pd.DataFrame(dict(d=d, d0=R0)).dropna(); j['pd'] = j.d - j.d0; j = j.join(C[['lat', 'lon']])
    ci = block_ci(j, 'd', nb=1000); cp = block_ci(j, 'pd', nb=1000)
    rows.append(dict(variant=v, cells=len(j), med=j.d.median(), lo=ci[0], hi=ci[1], med_ref=j.d0.median(),
                     paired=j.pd.median(), p_lo=cp[0], p_hi=cp[1], shorter=(j.d < 0).mean()))
T = pd.DataFrame(rows).round(2); T.to_csv(f'{D}/res/fixed_{ref}_cap{cap}' + ('' if YEND == 2022 else f'_to{YEND}') + ('' if PAT == '*' else '_' + PAT.replace('*', '')) + '.csv', index=False); print(T.to_string())
