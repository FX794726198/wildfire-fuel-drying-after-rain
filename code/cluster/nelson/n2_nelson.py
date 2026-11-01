# Run the Nelson (2000) dead-fuel-moisture model (NFDRS4 implementation, 1-h and 10-h sticks) hourly 1979-2022 at a chunk of cells.
# Output per cell and LOCAL day: stick moisture (g/g) at 14:00 local solar time and the minimum over the local day.
# Usage: python n2_nelson.py <hourly dir> <libdfm.so> <chunk index> <n chunks> <outdir>
import numpy as np, pandas as pd, ctypes, sys, os
hd, lib, k, K, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]; os.makedirs(out, exist_ok=True)
L = ctypes.CDLL(lib); f = np.ctypeslib.ndpointer(np.float32, flags='C_CONTIGUOUS')
L.run_dfm.argtypes = [ctypes.c_int, ctypes.c_int, f, f, f, f, f]
C = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cells_sel.csv')); idx = np.array_split(np.arange(len(C)), K)[k]
Y = range(1979, 2023)
parts = {v: [] for v in 'T RH SW R'.split()}
for y in Y:
    z = np.load(f'{hd}/h_{y}.npz')
    for v in parts: parts[v].append(z[v][idx, 24:-24])           # drop the 24-h margins -> contiguous UTC hours
H = {v: np.ascontiguousarray(np.concatenate(parts[v], 1)) for v in parts}
nh = H['T'].shape[1]; nd = nh // 24
res = {}
for lag in (1, 10):
    M = np.zeros((len(idx), nh), np.float32)
    for j in range(len(idx)):
        o = np.zeros(nh, np.float32); L.run_dfm(lag, nh, H['T'][j], H['RH'][j], H['SW'][j], H['R'][j], o); M[j] = o
    lon = C.lon.values[idx]; off = np.rint(-lon / 15).astype(int)                 # UTC hour of local midnight
    d = np.arange(1, nd - 1)
    h14 = d[None, :] * 24 + off[:, None] + 14
    res[f'm{lag}_14'] = np.take_along_axis(M, h14, 1)
    hh = (d[None, :] * 24 + off[:, None])[:, :, None] + np.arange(24)[None, None, :]
    res[f'm{lag}_min'] = np.take_along_axis(M, hh.reshape(len(idx), -1), 1).reshape(len(idx), len(d), 24).min(-1)
np.savez_compressed(f'{out}/nelson_{k:03d}.npz', cells=C.cell.values[idx], day0=1, **res)
print('chunk', k, len(idx), {kk: float(np.nanmean(v)) for kk, v in res.items()}, flush=True)
