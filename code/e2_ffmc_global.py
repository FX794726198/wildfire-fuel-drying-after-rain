# Global FFMC (Van Wagner 1987), vectorised over cells, and wet-fuel days after fire-season rains (same rules as the CEE paper).
# Inputs: daily arrays [cell, day] for 1979..Y1 from e1 (ERA5) or e3 (CMIP6). Output: per-cell period summary + per-event table.
import numpy as np, pandas as pd

def step(f, t, h, w, r):
    """One FFMC day for arrays; returns new FFMC (inputs assumed finite)."""
    mo = 147.2 * (101 - f) / (59.5 + f)
    rf = np.maximum(r - 0.5, 1e-9); wet = r > 0.5
    add = 42.5 * rf * np.exp(-100 / (251 - mo)) * (1 - np.exp(-6.93 / rf))
    add = add + np.where(mo > 150, 0.0015 * np.maximum(mo - 150, 0) ** 2 * np.sqrt(rf), 0)
    mo = np.where(wet, np.minimum(mo + add, 250), mo)
    ed = 0.942 * h ** 0.679 + 11 * np.exp((h - 100) / 10) + 0.18 * (21.1 - t) * (1 - np.exp(-0.115 * h))
    ew = 0.618 * h ** 0.753 + 10 * np.exp((h - 100) / 10) + 0.18 * (21.1 - t) * (1 - np.exp(-0.115 * h))
    ko = 0.424 * (1 - (h / 100) ** 1.7) + 0.0694 * np.sqrt(w) * (1 - (h / 100) ** 8)
    k1 = 0.424 * (1 - ((100 - h) / 100) ** 1.7) + 0.0694 * np.sqrt(w) * (1 - ((100 - h) / 100) ** 8)
    kd = ko * 0.581 * np.exp(0.0365 * t); kw = k1 * 0.581 * np.exp(0.0365 * t)
    m = np.where(mo > ed, ed + (mo - ed) * 10 ** (-kd), np.where(mo < ew, ew - (ew - mo) * 10 ** (-kw), mo))
    return 59.5 * (250 - m) / (147.2 + m)

def ffmc(T, H, W, P, F0=85.):
    """T degC, H %, W km/h, P mm/day; arrays [cell, day]. NaN input -> FFMC NaN, state carried."""
    n, nd = T.shape; F = np.full((n, nd), np.nan, np.float32); f = np.full(n, F0)
    for i in range(nd):
        t, h, w, r = T[:, i].astype(float), H[:, i].astype(float), W[:, i].astype(float), P[:, i].astype(float)
        ok = np.isfinite(t) & np.isfinite(h) & np.isfinite(r); w = np.where(np.isfinite(w), w, 12.)
        fn = step(f, np.where(ok, t, 20.), np.clip(np.where(ok, h, 50.), 1, 100), w, np.where(ok, r, 0.))
        f = np.where(ok, fn, f); F[:, i] = np.where(ok, fn, np.nan)
    return F

def rh(T, Td):
    es = lambda x: np.exp(17.625 * x / (243.04 + x)); return np.clip(100 * es(Td) / es(T), 1, 100)

def detrend(X, mon, yrs, ref):
    """Remove the per-cell, per-calendar-month linear trend in X (cells x days), referenced to year `ref`."""
    Y = X.copy()
    for m in range(1, 13):
        k = mon == m; x = yrs[k] - ref; xm = x.mean()
        v = X[:, k]; vm = np.nanmean(v, 1, keepdims=True)
        b = np.nansum((v - vm) * (x - xm), 1) / np.nansum(np.where(np.isfinite(v), (x - xm) ** 2, 0), 1)
        Y[:, k] = v - b[:, None] * x[None, :]
    return Y

def events(F, P, mon, yr, thr=10., rt=85., season_thr=80., ref=None, smon=None, pre=3):
    """Event table and fire-season months (FFMC climatology >= season_thr; over years in `ref` (y0, y1) if given).
    smon: optional fixed fire-season months [cell, 12]; pre: days before the rain that must have valid FFMC."""
    n, nd = F.shape; kr = np.ones(nd, bool) if ref is None else (yr >= ref[0]) & (yr <= ref[1])
    if smon is None:
        mF = np.stack([np.nanmean(F[:, (mon == m) & kr], 1) for m in range(1, 13)], 1)      # cell x 12
        smon = mF >= season_thr
    season = smon[:, mon - 1]
    ok = np.isfinite(F)
    rows = []
    for c in range(n):
        if season[c].sum() < 1000: continue
        cand = np.where(season[c] & (P[c] >= thr))[0]
        cand = cand[(cand >= pre) & (cand + 31 < nd)]
        for i in cand:
            if not ok[c, i - pre:i + 1].all() or not F[c, i - 1] >= rt: continue
            seg = F[c, i + 1:i + 31]; kk = np.where(seg >= rt)[0]; rec = kk[0] + 1 if len(kk) else 31
            if np.isnan(seg[:min(rec, 30)]).any(): continue
            fut = P[c, i + 1:i + 11]
            rows.append((c, yr[i], mon[i], i, rec, P[c, i], F[c, i - 1], F[c, i], np.nansum(fut >= 1) > 0, np.nanmax(fut) >= thr))
    ev = pd.DataFrame(rows, columns=['cell', 'y', 'm', 'i', 'rec', 'p', 'fpre', 'f0', 'further1', 'further10'])
    return ev, smon

def summarise(ev, F, P, H, T, smon, mon, yr, split, y0, y1):
    out = []
    for c, e in ev.groupby('cell'):
        E_, L_ = e[e.y < split], e[e.y >= split]
        if len(e) < 30 or len(E_) < 10 or len(L_) < 10: continue
        s = smon[c, mon - 1]; se, sl = s & (yr < split), s & (yr >= split)
        ok = np.isfinite(F[c])
        nye, nyl = split - y0, y1 - split + 1
        out.append(dict(cell=c, ne=len(E_), nl=len(L_), rec_e=E_.rec.mean(), rec_l=L_.rec.mean(),
                        rec15_e=E_.rec.clip(upper=15).mean(), rec15_l=L_.rec.clip(upper=15).mean(),
                        p_e=E_.p.median(), p_l=L_.p.median(), nev_e=len(E_) / nye, nev_l=len(L_) / nyl,
                        fpre_e=E_.fpre.median(), fpre_l=L_.fpre.median(), f0_e=E_.f0.median(), f0_l=L_.f0.median(),
                        H_e=np.nanmean(H[c, se]), H_l=np.nanmean(H[c, sl]), T_e=np.nanmean(T[c, se]), T_l=np.nanmean(T[c, sl]),
                        wet_e=(F[c, se & ok] < 85).sum() / nye, wet_l=(F[c, sl & ok] < 85).sum() / nyl,
                        rd_e=np.nansum(P[c, se] >= 1) / nye, rd_l=np.nansum(P[c, sl] >= 1) / nyl,
                        nseas=int(smon[c].sum()), months=''.join('1' if v else '0' for v in smon[c])))
    R = pd.DataFrame(out); R['d'] = R.rec_l / R.rec_e - 1; return R

def norain_rec(ev, F, T, H, W, rt=85.):
    """Recovery of each event with the FFMC run forward from the rain day under observed T/RH/wind but no further rain."""
    c, i = ev.cell.values, ev.i.values; f = F[c, i].astype(float); rec = np.full(len(ev), 31)
    done = np.zeros(len(ev), bool)
    for k in range(1, 31):
        t, h, w = T[c, i + k].astype(float), H[c, i + k].astype(float), W[c, i + k].astype(float)
        good = np.isfinite(t) & np.isfinite(h)
        f = np.where(good, step(f, np.where(good, t, 20.), np.clip(np.where(good, h, 50.), 1, 100), np.where(np.isfinite(w), w, 12.), 0.), f)
        hit = (~done) & (f >= rt); rec[hit] = k; done |= hit
    return rec
