# Global ERA5 (WeatherBench2 hourly 1.5 deg, conservative regrid) -> daily local-noon fire weather at land cells, 1979-2022.
# Per local day D and cell: T, Td, wind at the UTC hour closest to 12:30 local solar time; rain = 24 h ending at that hour
# (FWI convention); also Tmax and minimum RH over the local calendar day (for comparison with CMIP6 tasmax/hursmin).
# Usage: python e1_era5_extract.py <year> <outdir>     (one process per year; ~3 GB read per year)
import numpy as np, requests, numcodecs, sys, os, time, datetime as dt
B = 'https://storage.googleapis.com/weatherbench2/datasets/era5/1959-2023_01_10-1h-240x121_equiangular_with_poles_conservative.zarr'
bl = numcodecs.Blosc(); S = requests.Session()
S.mount('https://', requests.adapters.HTTPAdapter(pool_maxsize=32))
def get(p):
    for k in range(8):
        try:
            r = S.get(f'{B}/{p}', timeout=120)
            if r.status_code == 200: return r.content
        except Exception: pass
        time.sleep(2 ** k)
    raise RuntimeError(p)
def arr(v, c): return np.frombuffer(bl.decode(get(f'{v}/{c}.0.0')), '<f4').reshape(48, 240, 121)
def cells(outdir):
    f = f'{outdir}/cells.npz'
    if os.path.exists(f): return np.load(f)
    lsm = np.frombuffer(bl.decode(get('land_sea_mask/0.0')), '<f4').reshape(240, 121)
    lat = np.frombuffer(bl.decode(get('latitude/0')), '<f8'); lon = np.frombuffer(bl.decode(get('longitude/0')), '<f8')
    LO, LA = np.meshgrid(lon, lat, indexing='ij')
    k = (lsm >= 0.5) & (LA >= -56) & (LA <= 72)
    io, ia = np.where(k); lo = np.where(LO[k] > 180, LO[k] - 360, LO[k])
    np.savez(f, io=io, ia=ia, lat=LA[k], lon=lo, lsm=lsm[k]); return np.load(f)
if __name__ == '__main__':
    y, out = int(sys.argv[1]), sys.argv[2]; os.makedirs(out, exist_ok=True)
    C = cells(out); io, ia, lon = C['io'], C['ia'], C['lon']; n = len(io)
    h0 = int((dt.datetime(y, 1, 1) - dt.datetime(1959, 1, 1)).total_seconds() // 3600)
    nd = 366 if y % 4 == 0 else 365
    off = np.rint(12.5 - lon / 15).astype(int)          # UTC hour of local noon relative to 00 UTC of the local day
    dstart = -np.rint(lon / 15).astype(int)             # UTC hour of local midnight
    a, b = h0 - 48, h0 + nd * 24 + 48                   # hour window to load
    ca, cb = a // 48, (b - 1) // 48
    V = {v: np.empty((n, (cb - ca + 1) * 48), np.float32) for v in ['2m_temperature', '2m_dewpoint_temperature', 'total_precipitation', '10m_u_component_of_wind', '10m_v_component_of_wind']}
    from concurrent.futures import ThreadPoolExecutor
    def fill(job):
        j, c, v = job; V[v][:, j * 48:(j + 1) * 48] = arr(v, c)[:, io, ia].T
    with ThreadPoolExecutor(int(os.environ.get('NTHR', 16))) as ex:
        list(ex.map(fill, [(j, c, v) for j, c in enumerate(range(ca, cb + 1)) for v in V]))
    base = ca * 48
    T, Td = V['2m_temperature'] - 273.15, V['2m_dewpoint_temperature'] - 273.15
    W = np.sqrt(V['10m_u_component_of_wind'] ** 2 + V['10m_v_component_of_wind'] ** 2) * 3.6
    P = np.maximum(V['total_precipitation'], 0) * 1000.; Pc = np.cumsum(P, axis=1)
    es = lambda x: np.exp(17.625 * x / (243.04 + x)); H = np.clip(100 * es(Td) / es(T), 1, 100)
    rows = np.arange(n)[:, None]; days = np.arange(nd)[None, :]
    noon = h0 - base + days * 24 + off[:, None]          # column index of the noon hour
    o = dict(T12=T[rows, noon], Td12=Td[rows, noon], W12=W[rows, noon], P24=Pc[rows, noon] - Pc[rows, noon - 24])
    ds = h0 - base + days * 24 + dstart[:, None]
    hh = ds[:, :, None] + np.arange(24)[None, None, :]
    o['Tmax'] = T[rows[:, :, None], hh].max(-1); o['RHmin'] = H[rows[:, :, None], hh].min(-1)
    np.savez_compressed(f'{out}/era5_{y}.npz', **{k: v.astype(np.float32) for k, v in o.items()})
    print(y, n, {k: float(np.nanmean(v)) for k, v in o.items()})
