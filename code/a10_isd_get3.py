# Download isd-lite files for candidate stations and keep, for each local calendar day, the maximum temperature and minimum relative
# humidity over the 3-hourly synoptic observations (00, 03, ..., 21 UTC) and their number, for comparison with CMIP6 tasmax/hursmin.
# Local day: UTC hours [24 d - round(lon/15), +24), as for ERA5 Tmax/RHmin (e1_era5_extract.py). Output: isd/st3/{sid}.npz.
import numpy as np, pandas as pd, requests, gzip, io, os, sys
from concurrent.futures import ThreadPoolExecutor
D = './isd'; S = pd.read_csv(f'{D}/stations_cand.csv', dtype={'USAF': str, 'WBAN': str})
t0 = pd.Timestamp('1979-01-01'); nd = len(pd.date_range('1979-01-01', '2022-12-31'))
L = {y: set(pd.read_csv(f'{D}/lists/{y}.tsv', sep='\t', names=['k', 'size'], dtype={'k': str}).k) for y in range(1979, 2023)}
Ses = requests.Session(); Ses.mount('https://', requests.adapters.HTTPAdapter(pool_maxsize=40))
es = lambda x: np.exp(17.625 * x / (243.04 + x))
def one(r):
    sid = f'{r.USAF}-{r.WBAN}'; f = f'{D}/st3/{sid}.npz'
    if os.path.exists(f): return
    dstart = -int(np.rint(r.lon / 15)); Tx = np.full(nd, -np.inf); Hn = np.full(nd, np.inf); n3 = np.zeros(nd, np.int16)
    for y in range(1978, 2023):                          # 1978 for local days that start before 00 UTC on 1 Jan 1979
        if f'{sid}-{y}' not in L.get(y, ()) and y >= 1979: continue
        q = None
        for k in range(5):
            try:
                q = Ses.get(f'https://noaa-isd-pds.s3.amazonaws.com/isd-lite/data/{y}/{sid}-{y}.gz', timeout=60)
                if q.status_code == 200: break
            except Exception: pass
        if q is None or q.status_code != 200: continue
        a = np.loadtxt(io.StringIO(gzip.decompress(q.content).decode()), dtype=np.int32, ndmin=2)
        if not len(a): continue
        a = a[(a[:, 3] % 3 == 0) & (a[:, 4] != -9999) & (a[:, 5] != -9999)]
        T, Td = a[:, 4] / 10., a[:, 5] / 10.
        g = (T > -60) & (T < 60) & (Td > -70) & (Td <= T + 0.5); a, T, Td = a[g], T[g], np.minimum(Td[g], T[g])
        hrs = ((pd.to_datetime(dict(year=a[:, 0], month=a[:, 1], day=a[:, 2])) - t0).dt.days.values * 24 + a[:, 3]) - dstart
        di = np.floor_divide(hrs, 24); k = (di >= 0) & (di < nd); di, T, H = di[k], T[k], np.clip(100 * es(Td[k]) / es(T[k]), 1, 100)
        np.maximum.at(Tx, di, T); np.minimum.at(Hn, di, H); np.add.at(n3, di, 1)
    Tx[~np.isfinite(Tx)] = np.nan; Hn[~np.isfinite(Hn)] = np.nan
    np.savez_compressed(f, Tx=Tx.astype(np.float32), Hn=Hn.astype(np.float32), n3=n3, dstart=dstart)
    print(sid, int((n3 >= 7).sum()), flush=True)
with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 24) as ex: list(ex.map(one, [r for r in S.itertuples()]))
