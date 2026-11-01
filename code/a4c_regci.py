# Cell-bootstrap 95% intervals of the within-cell dose-response slope for each AR6 region (Fig. S3b).
import numpy as np, pandas as pd, sys
tag = sys.argv[1] if len(sys.argv) > 1 else 'modis'
B = pd.read_csv(f'res/fire_events_{tag}.csv.gz'); B['rc'] = B.rec.clip(upper=15)
def slope(b):
    x = b.rc - b.groupby('cell').rc.transform('mean'); y = b.avoid - b.groupby('cell').avoid.transform('mean')
    return float((x * y).sum() / (x * x).sum())
rng = np.random.default_rng(3); rows = []
for r, b in B.groupby('ar6'):
    if len(b) < 300 or b.cell.nunique() <= 5: continue
    G = {c: g for c, g in b.groupby('cell')}; cl = list(G)
    bs = [slope(pd.concat([G[c] for c in rng.choice(cl, len(cl))])) for _ in range(300)]
    lo, hi = np.percentile(bs, [2.5, 97.5]); rows.append(dict(ar6=r, n=len(b), cells=len(cl), slope=slope(b), lo=lo, hi=hi))
pd.DataFrame(rows).to_csv(f'res/fire_dose_{tag}_regions_ci.csv', index=False); print(pd.DataFrame(rows).round(3).to_string())
