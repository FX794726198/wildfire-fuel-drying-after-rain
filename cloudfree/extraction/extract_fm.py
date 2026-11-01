# Daily MODIS fire-mask counts (MOD14A1/MYD14A1 C6.1, 1 km) per ERA5 1.5-deg cell of the cloud-free analysis, with observation state.
# For every cell, day and sensor: pixel counts by QA land/water state (bits 0-1: 0 water, 1 coast, 2 land, 3 unused) x FireMask
# class (0-9: 0-2 not processed, 3 water, 4 cloud, 5 clear land without fire, 6 unknown, 7-9 fire low/nominal/high) = 40 counts.
# One output per tile and sensor (skipped if present), so the job can be restarted.
# Usage: python extract_fm.py urls.csv cells_cloudfree.csv cells_all.csv OUTDIR NPROC [MAXFILES]
import sys, os, subprocess, time, numpy as np, pandas as pd
from multiprocessing import Pool
URLS, CP, CA, OUT, NP = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
MAXF = int(sys.argv[6]) if len(sys.argv) > 6 else 0
TMP = os.environ.get('TMPDIR', '/tmp'); os.makedirs(OUT, exist_ok=True)
D0 = pd.Timestamp('2002-12-01'); ND = (pd.Timestamp('2022-12-31') - D0).days + 1
CPL = pd.read_csv(CP); CAL = pd.read_csv(CA); NC = len(CPL)
lut = np.full((121, 240), -1, np.int32); pidx = {c: k for k, c in enumerate(CPL.cell)}
for c, la, lo in zip(CAL.cell, CAL.lat, CAL.lon):
    if c in pidx: lut[int(round((la + 90) / 1.5)), int(round((lo % 360) / 1.5)) % 240] = pidx[c]
R = 6371007.181; T = 1111950.5197665; PS = T / 1200

def cellmap(h, v):
    """Cell index (0..NC-1, -1 outside) of every pixel of sinusoidal tile h, v (pixel centres)."""
    i = np.arange(1200) + 0.5; x = (h - 18) * T + i * PS; y = (9 - v) * T - i * PS
    X, Y = np.meshgrid(x, y); lat = np.degrees(Y / R)
    with np.errstate(divide='ignore', invalid='ignore'): lon = np.degrees(X / (R * np.cos(np.radians(lat))))
    m = np.full(X.shape, -1, np.int32); ok = np.isfinite(lon) & (np.abs(lon) <= 180)
    m[ok] = lut[np.rint((lat[ok] + 90) / 1.5).astype(int).clip(0, 120), np.rint((lon[ok] % 360) / 1.5).astype(int) % 240]
    return m

def read(f):
    try:
        from pyhdf.SD import SD, SDC
        s = SD(f, SDC.READ); fm = s.select('FireMask')[:]; qa = s.select('QA')[:]; s.end(); return fm, qa
    except ImportError:
        from osgeo import gdal
        g = lambda sds: gdal.Open(f'HDF4_EOS:EOS_GRID:"{f}":MODIS_Grid_Daily_Fire:{sds}').ReadAsArray()
        return g('FireMask'), g('QA')

CM = {}
def work(row):
    sensor, h, v, name, url = row
    if (h, v) not in CM: CM[(h, v)] = cellmap(h, v)
    m = CM[(h, v)]; sel = m >= 0
    if not sel.any(): return name, None, 'nocell'
    pid = os.getpid(); f = f'{TMP}/fm_{pid}.hdf'; ck = f'{TMP}/ck_{pid}.txt'
    for k in range(6):
        rc = subprocess.run(['curl', '-sS', '-n', '-L', '--fail', '-c', ck, '-b', ck, '--max-time', '600', '-o', f, url], capture_output=True)
        if rc.returncode == 0: break
        time.sleep(5 * (k + 1))
    else:
        return name, None, 'fail ' + rc.stderr.decode()[:200]
    try:
        fm, qa = read(f)
    except Exception as e:
        return name, None, 'read ' + str(e)[:200]
    if fm.ndim == 2: fm, qa = fm[None], qa[None]
    a = name.split('.')[1]; start = pd.Timestamp(int(a[1:5]), 1, 1) + pd.Timedelta(days=int(a[5:8]) - 1)
    mi = m[sel]; res = []
    for k in range(fm.shape[0]):
        di = (start - D0).days + k
        if di < 0 or di >= ND: continue
        if (start + pd.Timedelta(days=k)).year != start.year: continue        # last granule of a year: no spill-over layers
        cls = fm[k][sel].astype(np.int64); ls = (qa[k][sel] & 3).astype(np.int64)
        ok = (cls >= 0) & (cls <= 9)
        cnt = np.bincount(mi[ok] * 40 + ls[ok] * 10 + cls[ok], minlength=NC * 40).reshape(NC, 40)
        res.append((di, cnt.astype(np.uint16)))
    return name, res, 'ok'

if __name__ == '__main__':
    U = pd.read_csv(URLS)
    for (sensor, h, v), g in U.groupby(['sensor', 'h', 'v']):
        fn = f'{OUT}/fm_{sensor}_h{h:02d}v{v:02d}.npz'
        if os.path.exists(fn): continue
        rows = list(g.itertuples(index=False, name=None))
        if MAXF: rows = rows[:MAXF]
        A = np.zeros((NC, ND, 40), np.uint16); got = np.zeros(ND, bool); log = []; t0 = time.time()
        with Pool(NP) as p:
            for name, res, st in p.imap_unordered(work, rows, chunksize=2):
                log.append((name, st))
                if res is None: continue
                for di, cnt in res:
                    A[:, di] += cnt; got[di] = True
        keep = np.where(A.reshape(NC, -1).any(1))[0]
        np.savez_compressed(fn + '.tmp.npz', counts=A[keep], cells=CPL.cell.values[keep], got=got, d0=str(D0.date()))
        os.replace(fn + '.tmp.npz', fn)
        pd.DataFrame(log, columns=['name', 'status']).to_csv(fn.replace('.npz', '_log.csv'), index=False)
        bad = sum(s != 'ok' for _, s in log)
        print(sensor, h, v, len(rows), 'files', bad, 'not ok', len(keep), 'cells', round(time.time() - t0), 's', flush=True)
