# Hourly ERA5 (WeatherBench2 1.5 deg) at the 1,704 analysis cells for one year, for the Nelson dead-fuel-moisture model.
# Saves T (degC), RH (0-1), SWdown (W m-2, net shortwave / 0.8), rain (cm per hour) as float32 [cell, hour], hours of the UTC year
# plus 24 h on each side. Usage: python n1_hourly.py <year> <outdir>
import numpy as np, requests, numcodecs, sys, os, time, datetime as dt, pandas as pd
from concurrent.futures import ThreadPoolExecutor
B = 'https://storage.googleapis.com/weatherbench2/datasets/era5/1959-2023_01_10-1h-240x121_equiangular_with_poles_conservative.zarr'
bl = numcodecs.Blosc(); S = requests.Session(); S.mount('https://', requests.adapters.HTTPAdapter(pool_maxsize=32))
def get(p):
    for k in range(8):
        try:
            r = S.get(f'{B}/{p}', timeout=120)
            if r.status_code == 200: return r.content
        except Exception: pass
        time.sleep(2 ** k)
    raise RuntimeError(p)
y, out = int(sys.argv[1]), sys.argv[2]; os.makedirs(out, exist_ok=True)
C = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cells_sel.csv')); io, ia = C.io.values, C.ia.values; n = len(C)
h0 = int((dt.datetime(y, 1, 1) - dt.datetime(1959, 1, 1)).total_seconds() // 3600); nh = (366 if y % 4 == 0 else 365) * 24
a, b = h0 - 24, h0 + nh + 24; ca, cb = a // 48, (b - 1) // 48
VN = ['2m_temperature', '2m_dewpoint_temperature', 'total_precipitation', 'mean_surface_net_short_wave_radiation_flux']
V = {v: np.empty((n, (cb - ca + 1) * 48), np.float32) for v in VN}
def fill(job):
    j, c, v = job; V[v][:, j * 48:(j + 1) * 48] = np.frombuffer(bl.decode(get(f'{v}/{c}.0.0')), '<f4').reshape(48, 240, 121)[:, io, ia].T
with ThreadPoolExecutor(12) as ex: list(ex.map(fill, [(j, c, v) for j, c in enumerate(range(ca, cb + 1)) for v in VN]))
s = slice(a - ca * 48, b - ca * 48)
T = V['2m_temperature'][:, s] - 273.15; Td = V['2m_dewpoint_temperature'][:, s] - 273.15
es = lambda x: np.exp(17.625 * x / (243.04 + x)); RH = np.clip(es(Td) / es(T), 0.01, 1.0)
np.savez(f'{out}/h_{y}.npz', T=T, RH=RH.astype(np.float32), SW=np.maximum(V['mean_surface_net_short_wave_radiation_flux'][:, s] / 0.8, 0),
         R=np.maximum(V['total_precipitation'][:, s], 0) * 100., h0=a)
print(y, n, float(T.mean()), float(RH.mean()), flush=True)
