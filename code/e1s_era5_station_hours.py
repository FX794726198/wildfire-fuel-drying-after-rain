# ERA5 (WeatherBench2 1.5 deg) 2 m T and Td at each ISD station's fixed synoptic UTC hour (the 3-hourly hour closest to 12:30 local
# solar time), for the station's ERA5 cell and each local day. Pairs ERA5 and station humidity at the same time of day.
# Usage: python e1s_era5_station_hours.py <stations csv> <outdir> <y0> <y1>
import numpy as np, pandas as pd, sys, os, time, datetime as dt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e1_era5_extract import arr
from concurrent.futures import ThreadPoolExecutor
S = pd.read_csv(sys.argv[1]); out = sys.argv[2]; os.makedirs(out, exist_ok=True)
C = np.load('./era5g/cells.npz'); io, ia = C['io'][S.cell.values], C['ia'][S.cell.values]
hs = (3 * np.rint((12.5 - S.lon.values / 15) / 3)).astype(int)        # 0..24 h after 00 UTC of the local day
for y in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
    f = f'{out}/era5s_{y}.npz'
    if os.path.exists(f): continue
    t0 = time.time(); h0 = int((dt.datetime(y, 1, 1) - dt.datetime(1959, 1, 1)).total_seconds() // 3600)
    nd = 366 if y % 4 == 0 else 365; a, b = h0, h0 + nd * 24 + 25; ca, cb = a // 48, (b - 1) // 48
    V = {v: np.empty((len(S), (cb - ca + 1) * 48), np.float32) for v in ['2m_temperature', '2m_dewpoint_temperature']}
    def fill(job):
        j, c, v = job; V[v][:, j * 48:(j + 1) * 48] = arr(v, c)[:, io, ia].T
    with ThreadPoolExecutor(int(os.environ.get('NTHR', 12))) as ex:
        list(ex.map(fill, [(j, c, v) for j, c in enumerate(range(ca, cb + 1)) for v in V]))
    col = h0 - ca * 48 + np.arange(nd)[None, :] * 24 + hs[:, None]; rows = np.arange(len(S))[:, None]
    np.savez_compressed(f, T=V['2m_temperature'][rows, col] - 273.15, Td=V['2m_dewpoint_temperature'][rows, col] - 273.15, hs=hs)
    print(y, round(time.time() - t0), flush=True)
