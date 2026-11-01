# ERA5 (WeatherBench2 hourly 1.5 deg) daily-mean 2 m temperature and relative humidity over the local calendar day at the 1,704
# analysis cells, for comparing observed fire-season changes with DAMIP monthly-mean (all-hour) model responses.
# Usage: python e1m_era5_dailymean.py <outdir> <y0> <y1>
import numpy as np, pandas as pd, sys, os, time, datetime as dt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e1_era5_extract import arr
from concurrent.futures import ThreadPoolExecutor
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
SEL = pd.read_csv('./code/cells_sel.csv').cell.values
C = np.load('./era5g/cells.npz'); io, ia, lon = C['io'][SEL], C['ia'][SEL], C['lon'][SEL]
dstart = -np.rint(lon / 15).astype(int); es = lambda x: np.exp(17.625 * x / (243.04 + x))
for y in range(int(sys.argv[2]), int(sys.argv[3]) + 1):
    f = f'{out}/era5m_{y}.npz'
    if os.path.exists(f): continue
    t0 = time.time(); h0 = int((dt.datetime(y, 1, 1) - dt.datetime(1959, 1, 1)).total_seconds() // 3600)
    nd = 366 if y % 4 == 0 else 365; a, b = h0 - 24, h0 + nd * 24 + 24; ca, cb = a // 48, (b - 1) // 48
    V = {v: np.empty((len(SEL), (cb - ca + 1) * 48), np.float32) for v in ['2m_temperature', '2m_dewpoint_temperature']}
    def fill(job):
        j, c, v = job; V[v][:, j * 48:(j + 1) * 48] = arr(v, c)[:, io, ia].T
    with ThreadPoolExecutor(int(os.environ.get('NTHR', 12))) as ex:
        list(ex.map(fill, [(j, c, v) for j, c in enumerate(range(ca, cb + 1)) for v in V]))
    col = h0 - ca * 48 + np.arange(nd)[None, :, None] * 24 + dstart[:, None, None] + np.arange(24)[None, None, :]
    rows = np.arange(len(SEL))[:, None, None]
    T = V['2m_temperature'][rows, col] - 273.15; Td = V['2m_dewpoint_temperature'][rows, col] - 273.15
    H = np.clip(100 * es(Td) / es(T), 1, 100)
    np.savez_compressed(f, Tm=T.mean(-1).astype(np.float32), Hm=H.mean(-1).astype(np.float32))
    print(y, round(time.time() - t0), flush=True)
