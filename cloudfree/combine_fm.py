# Combine per tile/sensor fire-mask counts (extract_fm.py) into arrays [cell, day, sensor] for fires (classes 8-9 and 7-9),
# valid observations (5, 7, 8, 9) and cloud (4), land and coast pixels only (QA bits 0-1 = 1 or 2).
# Usage: python combine_fm.py <dir with fm_*.npz> <cells_cloudfree.csv> <out fm.npz>
import numpy as np, pandas as pd, glob, sys, os
D, CP, OUT = sys.argv[1:4]; pc = pd.read_csv(CP).cell.values; pidx = {c: k for k, c in enumerate(pc)}
f89 = f79 = valid = cloud = None; got = {}
for f in sorted(glob.glob(f'{D}/fm_*.npz')):
    s = 0 if 'MOD14A1' in f else 1; z = np.load(f); A = z['counts'].astype(np.float32); nd = A.shape[1]
    if f89 is None: f89, f79, valid, cloud = (np.zeros((len(pc), nd, 2), np.float32) for _ in range(4))
    L = A.reshape(A.shape[0], nd, 4, 10)[:, :, 1:3, :].sum(2)          # land + coast
    k = np.array([pidx[c] for c in z['cells']])
    f89[k, :, s] += L[..., 8] + L[..., 9]; f79[k, :, s] += L[..., 7] + L[..., 8] + L[..., 9]
    valid[k, :, s] += L[..., 5] + L[..., 7] + L[..., 8] + L[..., 9]; cloud[k, :, s] += L[..., 4]
    got[os.path.basename(f)] = int(z['got'].sum())
np.savez_compressed(OUT, f89=f89, f79=f79, valid=valid, cloud=cloud)
print(len(got), 'files', 'days with data per file: min', min(got.values()), 'max', max(got.values()))
