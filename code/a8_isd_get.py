# Download isd-lite files for candidate stations and keep T/Td at each station's fixed synoptic UTC hour for each local day.
# Output: isd/station_{sid}.npz with T, Td [days 1979-2022] (NaN where missing).
import numpy as np, pandas as pd, requests, gzip, io, os, sys
from concurrent.futures import ThreadPoolExecutor
D = './isd'; S = pd.read_csv(f'{D}/stations_cand.csv', dtype={'USAF': str, 'WBAN': str})
dates = pd.date_range('1979-01-01', '2022-12-31'); nd = len(dates); os.makedirs(f'{D}/st', exist_ok=True)
L = {}
for y in range(1979, 2023):
    t = pd.read_csv(f'{D}/lists/{y}.tsv', sep='\t', names=['k', 'size'], dtype={'k': str}); L[y] = set(t.k)
Ses = requests.Session(); Ses.mount('https://', requests.adapters.HTTPAdapter(pool_maxsize=40))
def one(r):
    sid = f'{r.USAF}-{r.WBAN}'; f = f'{D}/st/{sid}.npz'
    if os.path.exists(f): return
    hs = int(3 * np.rint((12.5 - r.lon / 15) / 3)); T = np.full(nd, np.nan, np.float32); Td = T.copy(); n = 0
    for y in range(1979, 2023):
        if f'{sid}-{y}' not in L[y]: continue
        for k in range(5):
            try:
                q = Ses.get(f'https://noaa-isd-pds.s3.amazonaws.com/isd-lite/data/{y}/{sid}-{y}.gz', timeout=60)
                if q.status_code == 200: break
            except Exception: pass
        else: continue
        a = np.loadtxt(io.StringIO(gzip.decompress(q.content).decode()), dtype=np.int32, ndmin=2)
        if not len(a): continue
        tt = pd.to_datetime(dict(year=a[:, 0], month=a[:, 1], day=a[:, 2])) + pd.to_timedelta(a[:, 3] - hs, 'h')
        ok = (tt.dt.hour.values == 0) & (a[:, 4] != -9999) & (a[:, 5] != -9999)       # exactly the station's synoptic hour
        di = ((tt[ok] - dates[0]).dt.days).values; g = (di >= 0) & (di < nd)
        T[di[g]] = a[ok, 4][g] / 10.; Td[di[g]] = a[ok, 5][g] / 10.; n += g.sum()
    np.savez_compressed(f, T=T, Td=Td, hs=hs)
    print(sid, n, flush=True)
with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 24) as ex: list(ex.map(one, [r for r in S.itertuples()]))
