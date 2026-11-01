# ERA5 (WeatherBench2 1.5 deg) daily maximum temperature and minimum relative humidity at each ISD candidate station's cell, sampled
# at the 3-hourly synoptic UTC hours of the station's local calendar day (the same hours as a10_isd_get3.py), for station comparison.
# Usage: python e1t_era5_txrn3.py <stations csv> <outdir> <y0> <y1>
import numpy as np, pandas as pd, sys, os, time, datetime as dt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e1_era5_extract import arr
from concurrent.futures import ThreadPoolExecutor
S = pd.read_csv(sys.argv[1]); out = sys.argv[2]; os.makedirs(out, exist_ok=True)
C = np.load('./era5g/cells.npz'); io, ia = C['io'][S.cell.values], C['ia'][S.cell.values]
dstart = -np.rint(S.lon.values / 15).astype(int)
hh = np.array([[h for h in range(24) if (ds + h) % 3 == 0] for ds in dstart])          # station x 8 hours of the local day
es = lambda x: np.exp(17.625 * x / (243.04 + x))
for y in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
    f = f'{out}/era5t_{y}.npz'
    if os.path.exists(f): continue
    t0 = time.time(); h0 = int((dt.datetime(y, 1, 1) - dt.datetime(1959, 1, 1)).total_seconds() // 3600)
    nd = 366 if y % 4 == 0 else 365; a, b = h0 - 24, h0 + nd * 24 + 24; ca, cb = a // 48, (b - 1) // 48
    V = {v: np.empty((len(S), (cb - ca + 1) * 48), np.float32) for v in ['2m_temperature', '2m_dewpoint_temperature']}
    def fill(job):
        j, c, v = job; V[v][:, j * 48:(j + 1) * 48] = arr(v, c)[:, io, ia].T
    with ThreadPoolExecutor(int(os.environ.get('NTHR', 12))) as ex:
        list(ex.map(fill, [(j, c, v) for j, c in enumerate(range(ca, cb + 1)) for v in V]))
    col = h0 - ca * 48 + np.arange(nd)[None, :, None] * 24 + dstart[:, None, None] + hh[:, None, :]   # station x day x 8
    rows = np.arange(len(S))[:, None, None]
    T = V['2m_temperature'][rows, col] - 273.15; Td = V['2m_dewpoint_temperature'][rows, col] - 273.15
    H = np.clip(100 * es(Td) / es(T), 1, 100)
    np.savez_compressed(f, Tx=T.max(-1).astype(np.float32), Hn=H.min(-1).astype(np.float32))
    print(y, round(time.time() - t0), flush=True)
