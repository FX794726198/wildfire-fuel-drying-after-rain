# Land-cover shares of each 1.5-degree cell of the cloud-free analysis from ESA WorldCover 2021 v200 (10 m map read at its 320 m
# overview, public on AWS). Usage: python wc_cover.py   (writes cell_cover.csv next to this script; set CURL_CA_BUNDLE if a proxy needs it)
import rasterio, os, numpy as np, pandas as pd, math, sys
HERE = os.path.dirname(os.path.abspath(__file__))
from rasterio.windows import from_bounds
CP = pd.read_csv(HERE + '/cells_cloudfree.csv')
def tname(la, lo): return f"{'N' if la >= 0 else 'S'}{abs(la):02d}{'E' if lo >= 0 else 'W'}{abs(lo):03d}"
cache = {}
def tile(la, lo):
    k = (la, lo)
    if k not in cache:
        u = f'/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{tname(la, lo)}_Map.tif'
        try:
            with rasterio.open(u) as r: cache[k] = r.read(1, out_shape=(1125, 1125))
        except Exception as e: cache[k] = None
    return cache[k]
rows = []
env = dict(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR')
if os.environ.get('CURL_CA_BUNDLE'): env['CURL_CA_BUNDLE'] = os.environ['CURL_CA_BUNDLE']
with rasterio.Env(**env):
    for _, c in CP.iterrows():
        s, n, w, e = c.lat - .75, c.lat + .75, c.lon - .75, c.lon + .75
        cnt = {}
        for la in range(int(math.floor(s / 3) * 3), int(math.ceil(n / 3) * 3), 3):
            for lo in range(int(math.floor(w / 3) * 3), int(math.ceil(e / 3) * 3), 3):
                a = tile(la, lo)
                if a is None: continue
                px = 3 / 1125                                    # degrees per pixel at this overview
                r0 = int(round((la + 3 - min(n, la + 3)) / px)); r1 = int(round((la + 3 - max(s, la)) / px))
                c0 = int(round((max(w, lo) - lo) / px)); c1 = int(round((min(e, lo + 3) - lo) / px))
                v, k = np.unique(a[r0:r1, c0:c1], return_counts=True)
                for vv, kk in zip(v, k): cnt[int(vv)] = cnt.get(int(vv), 0) + int(kk)
        land = sum(k for v, k in cnt.items() if v not in (0, 80))
        f = lambda v: cnt.get(v, 0) / land if land else np.nan
        rows.append(dict(cell=c.cell, ar6=c.ar6, tree=f(10), shrub=f(20), grass=f(30), crop=f(40), built=f(50), bare=f(60), wetland=f(90) + f(95), land_px=land))
        print(len(rows), end=' ', flush=True)
D = pd.DataFrame(rows); D.to_csv(HERE + '/cell_cover.csv', index=False)
print(); print(D.groupby('ar6')[['tree', 'shrub', 'grass', 'crop']].median().round(2)); print('tiles', len(cache), 'missing', sum(v is None for v in cache.values()))
