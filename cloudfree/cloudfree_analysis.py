# Fire activity on cloud-free land after soaking rains, by the reconstructed fuel drying time (Fig. 2e-g, Fig. S3c-e,g-l, Section S3).
# Five regions (SAM, SES, ESAF, NAU, CAF; 347 cells), 2003-2022. Inputs: daily MODIS fire-mask counts per cell (cloudfree/fm.npz, from
# combine_fm.py), ERA5 FFMC (out_era5/ffmc_obs.npz), ERA5 24-hour rain (era5g/era5_<year>.npz; cloudfree/era5g_p24 holds it for these
# cells), the rain events of the main analysis (out_era5/events_obs.csv.gz) and the tree cover of each cell (cell_cover.csv, wc_cover.py).
# Usage (from the package root): python cloudfree/cloudfree_analysis.py <group> [--boot 1000] [--only 'tag1|tag2'] [--name NAME]
#   groups: pooled, regions, periods, leave_one_out, conditions, cover, all (every result), check (fastglm against statsmodels),
#           merge (collects res/*.json into res/cloudfree_results.json and res/cloudfree_summary.csv)
# The Poisson fits use fastglm.py (numpy); with FASTGLM=0 they use statsmodels instead (same estimate, much slower).
import sys, os, json, glob, argparse, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE) + '/'
REG = ['SAM', 'SES', 'ESAF', 'NAU', 'CAF']
ap = argparse.ArgumentParser(); ap.add_argument('group'); ap.add_argument('--boot', type=int, default=1000)
ap.add_argument('--out', default=HERE + '/res'); ap.add_argument('--only', default=''); ap.add_argument('--name', default='')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

ORDER = ['five regions: main', 'five regions: clear>=0.7', 'five regions: ten-day baseline', 'five regions: confidence 7-9',
         'five regions: Terra', 'five regions: Aqua'] + [f'region {r}' for r in REG] + ['period 2003-2012', 'period 2013-2022'] + \
        [f'without {r}' for r in REG] + ['season: first two-thirds', 'season: last third', 'later rain: below 60 mm',
         'later rain: 60 mm or more', 'tree cover: below 26%', 'tree cover: 26% or more']
if a.group == 'merge':
    got = {}
    for f in sorted(glob.glob(a.out + '/*.json')):
        x = json.load(open(f))
        if os.path.basename(f) == 'cloudfree_results.json' or not isinstance(x, list): continue
        for r in x: got.setdefault(r['tag'], r)
    R = [got[t] for t in ORDER if t in got]
    json.dump(R, open(a.out + '/cloudfree_results.json', 'w'), indent=1)
    f2 = lambda x: np.round(x, 3)
    pd.DataFrame([dict(tag=r['tag'], events_1d=r['n_events'][0], events_2d=r['n_events'][1], events_3_15d=r['n_events'][2],
                       T50_1d=f2(r['T50'][0]), T50_2d=f2(r['T50'][1]), T50_3_15d=f2(r['T50'][2]), dT50=f2(r['dT50']),
                       dT50_lo=f2(r['dT50_ci'][0]), dT50_hi=f2(r['dT50_ci'][1]), rate_ratio_wet=f2(r['M3_rr_wet']),
                       rr_lo=f2(r['M3_ci'][0]), rr_hi=f2(r['M3_ci'][1]), level_days_8_15=f2(r['level_8_15'])) for r in R]
                 ).to_csv(a.out + '/cloudfree_summary.csv', index=False)
    print(len(R), 'results; missing:', [t for t in ORDER if t not in got]); sys.exit()
if os.environ.get('FASTGLM', '1') == '1':
    sys.path.insert(0, HERE); import fastglm; fastglm.install(sm)

# ---- inputs
CP = pd.read_csv(HERE + '/cells_cloudfree.csv'); pc = CP.cell.values
d79 = pd.date_range('1979-01-01', '2022-12-31'); D0 = pd.Timestamp('2002-12-01'); off = (D0 - d79[0]).days
dates = d79[off:]; ND = len(dates)
Z = np.load(HERE + '/fm.npz')                                       # each [cell, day from 2002-12-01, sensor 0 Terra / 1 Aqua]
F89, F79, V, CL = (Z[k] for k in ('f89', 'f79', 'valid', 'cloud'))  # fires of nominal/high and of all confidence, clear land, cloud
if os.path.exists(ROOT + 'out_era5/ffmc_obs.npz'):
    FF = np.load(ROOT + 'out_era5/ffmc_obs.npz'); FFMC = FF['F'][pc][:, off:].astype(np.float32); smon = FF['smon'][pc]
else:                                                               # subset for the 347 cells (cloudfree/ffmc_p347.npz)
    FF = np.load(HERE + '/ffmc_p347.npz'); ix = pd.Index(FF['cell']).get_indexer(pc); assert (ix >= 0).all()
    FFMC = FF['F'][ix][:, off:].astype(np.float32); smon = FF['smon'][ix]
P = np.concatenate([np.load(ROOT + f'era5g/era5_{y}.npz')['P24'][pc] for y in range(1979, 2023)], 1)[:, off:]

def combine(sens):
    """Fire and clear-land counts of one satellite or both; cell-days with < 20% of the observed land clear count as missing."""
    s = slice(None) if sens == 'both' else (0 if sens == 'terra' else 1)
    f89, f79, v, cl = (np.nansum(x[:, :, s], -1) if sens == 'both' else x[:, :, s] for x in (F89, F79, V, CL))
    with np.errstate(invalid='ignore', divide='ignore'): cf = v / (v + cl)
    miss = ~(cf >= 0.2) | ~(v > 0)
    return f89, f79, v, cf, miss

# ---- events: rains of the main analysis in these cells, 2003-2022, drying time <= 15 d, no further rain >= 10 mm within 15 days
E = pd.read_csv(ROOT + 'out_era5/events_obs.csv.gz'); E = E[E.cell.isin(pc) & (E.y >= 2003) & (E.rec <= 15)].copy()
E['k'] = E.cell.map({c: k for k, c in enumerate(pc)}); E['d'] = E.i - off
E = E[(E.d - 30 >= 0) & (E.d + 21 < ND)]
E['iso'] = [np.nanmax(P[k, d + 1:d + 16]) < 10 for k, d in zip(E.k, E.d)]
E = E[E.iso].merge(CP[['cell', 'lat', 'lon', 'ar6']], on='cell').reset_index(drop=True)
E['cls'] = np.select([E.rec <= 1, E.rec <= 2], [0, 1], 2)            # drying time 1 / 2 / 3-15 days
def third(k, m):
    """Position of the event month among the cell's fire-season months (circular order starting after the longest gap)."""
    s = np.where(smon[k])[0]
    if len(s) == 0: return 1
    start = next((j for j in range(12) if smon[k][j] and not smon[k][j - 1]), s[0])
    order = [(start + j) % 12 for j in range(12) if smon[k][(start + j) % 12]]
    pos = order.index(m - 1) if (m - 1) in order else 0
    return min(2, int(3 * pos / len(order)))
E['third'] = [third(k, m) for k, m in zip(E.k, E.m)]
E['rain16_45'] = [np.nansum(P[k, d + 16:d + 46]) for k, d in zip(E.k, E.d)]
E = E.merge(pd.read_csv(HERE + '/cell_cover.csv')[['cell', 'tree']], on='cell', how='left')
# dry-fuel days: FFMC >= 85 and no rain >= 10 mm in the previous 15 days
P10 = P >= 10; prior = np.zeros_like(P10)
for L in range(1, 16): prior[:, L:] |= P10[:, :-L]
DRY = (FFMC >= 85) & ~prior
T = np.arange(-10, 21)

def use(sub):
    """Select the events of one analysis and reset the resampling blocks (5 degrees x 1 year) and the random generator."""
    global EV, cl_id, cl_inv, rng
    EV = sub.reset_index(drop=True)
    blk = (np.floor(EV.lat / 5).astype(int) * 1000 + np.floor(EV.lon / 5).astype(int)) * 10000 + EV.y
    cl_id, cl_inv = np.unique(blk, return_inverse=True); rng = np.random.default_rng(20261028)
def event_arrays(f, v, miss, base):
    """Fire pixels and expected fire pixels (clear pixels x reference rate) on days -10..20 around each event."""
    n = len(EV); Fe = np.zeros((n, len(T))); VR = np.zeros((n, len(T))); R = np.full(n, np.nan)
    for e, (k, d) in enumerate(zip(EV.k, EV.d)):
        if base == 'season':                                          # dry-fuel days in the 20 days before the rain
            idx = np.arange(d - 20, d + 21); idx = idx[(idx >= 0) & (idx < ND) & ((idx < d) | (idx > d + 20))]
            idx = idx[DRY[k, idx] & ~miss[k, idx]]
            if len(idx) >= 10 and f[k, idx].sum() >= 5: R[e] = f[k, idx].sum() / v[k, idx].sum()
        else:                                                         # all observed days in the ten days before the rain
            idx = np.arange(d - 10, d); idx = idx[~miss[k, idx]]
            if len(idx) >= 5 and f[k, idx].sum() > 0: R[e] = f[k, idx].sum() / v[k, idx].sum()
        tt = d + T; ok = ~miss[k, tt]
        Fe[e] = np.where(ok, f[k, tt], 0); VR[e] = np.where(ok, v[k, tt], 0) * (R[e] if np.isfinite(R[e]) else 0)
    return Fe, VR, R
def strata(R, use_):
    """Weighting strata: region x third of the fire season x rain amount half x reference-rate half (halves by regional median)."""
    s = pd.DataFrame(dict(reg=EV.ar6, th=EV.third, cls=EV.cls))
    s['pq'] = (EV.p > EV.groupby('ar6').p.transform('median')).astype(int)
    Rs = pd.Series(R); s['rq'] = (Rs > Rs.groupby(EV.ar6).transform('median')).astype(int)
    s['st'] = s.reg.astype(str) + s.th.astype(str) + s.pq.astype(str) + s.rq.astype(str)
    s.loc[~use_, 'st'] = None
    return s
def weights(s, m):
    """Direct standardisation: each class's stratum distribution = the pooled one; m = event multiplicities (bootstrap).
    Events in strata that lack one of the classes get weight 0."""
    d = pd.DataFrame(dict(st=s.st, cls=s.cls, m=m)); d = d[d.st.notna()]
    tot = d.groupby('st').m.sum() / d.m.sum(); byc = d.groupby(['cls', 'st']).m.sum() / d.groupby('cls').m.sum()
    w = np.zeros(len(s)); ok = s.st.notna().values
    key = list(zip(s.cls[ok], s.st[ok])); w[ok] = [tot[st] / byc[(c, st)] if byc.get((c, st), 0) > 0 else 0 for c, st in key]
    present = d.groupby('st').cls.nunique(); full = s.st.map(present).fillna(0).values == 3
    return np.where(full, w, 0.0)
def curves(Fe, VR, w, m):
    out = []
    for c in range(3):
        sel = (EV.cls.values == c); ww = (w * m)[sel][:, None]
        num = (ww * Fe[sel]).sum(0); den = (ww * VR[sel]).sum(0)
        out.append(num / np.where(den > 0, den, np.nan))
    return np.array(out)                                               # fire activity relative to the reference [class, day]
def smooth(A): return np.array([np.convolve(np.nan_to_num(x, nan=np.nanmean(x)), np.ones(3) / 3, 'same') for x in A])
def tx(A, thr):
    """First time (linear interpolation) at or after day 1 when the three-day running mean reaches thr; 21 if never."""
    S = smooth(A); res = []
    for x in S:
        j = np.where((T >= 1) & (x >= thr))[0]
        if not len(j): res.append(21.0); continue
        j = j[0]
        if T[j] == 1 or x[j - 1] >= thr: res.append(float(T[j])); continue
        res.append(float(T[j - 1] + (thr - x[j - 1]) / (x[j] - x[j - 1])))
    return np.array(res, float)
def run(sens='both', conf='89', base='season', minclear=0.2, boot=a.boot, m3boot=a.boot, tag='main'):
    f89, f79, v, cf, miss = combine(sens); f = f89 if conf == '89' else f79
    if minclear > 0.2: miss = miss | ~(cf >= minclear)
    Fe, VR, R = event_arrays(f, v, miss, base)
    use_ = np.isfinite(R) & (R > 0); s = strata(np.where(use_, R, np.nan), use_)
    m1 = np.ones(len(EV)); w = weights(s, m1)
    A = curves(Fe, VR, w, m1); T50, T75 = tx(A, 0.5), tx(A, 0.75)
    n_cls = [int(((w > 0) & (EV.cls.values == c)).sum()) for c in range(3)]
    dropped = float(1 - (w > 0).sum() / max(use_.sum(), 1))
    CFc = []                                                          # cloud-free fraction around the events, same weights
    for c in range(3):
        sel = (EV.cls.values == c) & (w > 0); vals = []
        for j, t in enumerate(T):
            tt = EV.d.values[sel] + t; kk = EV.k.values[sel]; x = cf[kk, tt]; ww = w[sel]; ok = np.isfinite(x)
            vals.append(np.sum(ww[ok] * x[ok]) / max(np.sum(ww[ok]), 1e-9))
        CFc.append(vals)
    CFc = np.array(CFc)
    jj = T >= 1; xa, ya = CFc[:, jj].ravel(), A[:, jj].ravel(); okr = np.isfinite(ya)
    r2 = float(np.corrcoef(xa[okr], ya[okr])[0, 1] ** 2) if okr.sum() > 3 else np.nan
    def aligned(m):                                                   # mean fire activity on days 1-5 after the fuels had dried
        S2 = []
        for c in range(3):
            sel = (EV.cls.values == c) & (w > 0); num = den = 0.0
            for e in np.where(sel)[0]:
                for s_ in range(1, 6):
                    t = EV.rec.values[e] + s_
                    if t <= 20: j = t + 10; num += w[e] * m[e] * Fe[e, j]; den += w[e] * m[e] * VR[e, j]
            S2.append(num / den if den > 0 else np.nan)
        return np.array(S2)
    rainA = np.nanmean(A[:, (T >= 2) & (T <= 8)], 1); dryA = aligned(m1)
    out = dict(tag=tag, sens=sens, conf=conf, base=base, minclear=minclear, n_events=n_cls, dropped_frac=dropped,
               A=A.round(4).tolist(), T50=T50.tolist(), T75=T75.tolist(), dT50=float(T50[2] - T50[0]), dT75=float(T75[2] - T75[0]),
               A_day1=A[:, T == 1].ravel().round(3).tolist(), clear=CFc.round(3).tolist(), r2_clear=r2,
               M2_rain_2_8=rainA.round(3).tolist(), M2_dry_1_5=dryA.round(3).tolist(),
               M2_sd_rain=float(np.nanstd(rainA)), M2_sd_dry=float(np.nanstd(dryA)))
    if boot:                                                          # resampling blocks of 5 degrees and one year
        bT50 = []; bT75 = []; bA1 = []
        for b in range(boot):
            mc = rng.multinomial(len(cl_id), np.ones(len(cl_id)) / len(cl_id)); m = mc[cl_inv].astype(float)
            wb = weights(s, m); Ab = curves(Fe, VR, wb, m); t5, t7 = tx(Ab, 0.5), tx(Ab, 0.75)
            bT50.append(t5[2] - t5[0]); bT75.append(t7[2] - t7[0]); bA1.append(Ab[:, T == 1].ravel())
        out['dT50_ci'] = np.percentile(bT50, [2.5, 97.5]).tolist(); out['dT75_ci'] = np.percentile(bT75, [2.5, 97.5]).tolist()
        out['A_day1_hi'] = np.percentile(np.array(bA1), 97.5, 0).round(3).tolist()
    # Poisson regression of daily fire pixels on days 1-15: days before the FFMC first returned to 85 vs later days
    rows = []
    for e in np.where(w > 0)[0]:
        k, d = EV.k[e], EV.d[e]
        for t in range(1, 16):
            if VR[e, t + 10] <= 0: continue
            rows.append((e, t, int(t < EV.rec[e]), Fe[e, t + 10], np.log(VR[e, t + 10]), np.log(EV.p[e]), cf[k, d + t], f'{EV.ar6[e]}_{EV.m[e]}', cl_inv[e]))
    M = pd.DataFrame(rows, columns=['e', 't', 'wet', 'y', 'off', 'lp', 'cf', 'rm', 'cl'])
    X = pd.get_dummies(M[['t', 'rm']].astype(str), drop_first=True).astype(float)
    X['wet'] = M.wet; X['lp'] = M.lp; X['cf'] = M.cf; X = sm.add_constant(X)
    def fit(idx):
        r = sm.GLM(M.y.values[idx], X.values[idx], family=sm.families.Poisson(), offset=M.off.values[idx]).fit()
        return float(np.exp(r.params[list(X.columns).index('wet')]))
    out['M3_rr_wet'] = fit(np.arange(len(M))); out['M3_n'] = int(len(M))
    if boot:
        rr = []; groups = {c: np.where(M.cl.values == c)[0] for c in np.unique(M.cl)}; gk = list(groups)
        for b in range(min(boot, m3boot)):
            pick = rng.choice(len(gk), len(gk)); idx = np.concatenate([groups[gk[i]] for i in pick])
            try: rr.append(fit(idx))
            except Exception: pass
        out['M3_ci'] = np.percentile(rr, [2.5, 97.5]).tolist()
    return out

def level(A): return np.nanmean(np.asarray(A, float)[:, (T >= 8) & (T <= 15)], 1)
variants = [('main', {}), ('clear>=0.7', dict(minclear=0.7)), ('ten-day baseline', dict(base='pre')),
            ('confidence 7-9', dict(conf='79')), ('Terra', dict(sens='terra')), ('Aqua', dict(sens='aqua'))]
mr = float(E.rain16_45.median()); mt = float(E.drop_duplicates('cell').tree.median())
J = dict(pooled=[('five regions: ' + t, E, kw) for t, kw in variants],
         regions=[('region ' + r, E[E.ar6 == r], {}) for r in REG],
         periods=[('period 2003-2012', E[E.y <= 2012], {}), ('period 2013-2022', E[E.y >= 2013], {})],
         leave_one_out=[('without ' + r, E[E.ar6 != r], {}) for r in REG],
         conditions=[('season: first two-thirds', E[E.third < 2], {}), ('season: last third', E[E.third == 2], {}),
                     (f'later rain: below {mr:.0f} mm', E[E.rain16_45 < mr], {}), (f'later rain: {mr:.0f} mm or more', E[E.rain16_45 >= mr], {})],
         cover=[(f'tree cover: below {100 * mt:.0f}%', E[E.tree < mt], {}), (f'tree cover: {100 * mt:.0f}% or more', E[E.tree >= mt], {})])
J['all'] = sum(J.values(), [])
J['check'] = [('region SAM', E[E.ar6 == 'SAM'], dict(m3boot=30)), ('five regions: main', E, dict(m3boot=3))]
jobs = J[a.group]
if a.only: jobs = [j for j in jobs if j[0] in a.only.split('|')]
res = []
for tag, sub, kw in jobs:
    use(sub); r = run(tag=tag, **kw); lv = level(r['A'])
    r['level_8_15'] = float(np.nanmean(lv)); r['level_8_15_by_class'] = lv.round(3).tolist(); r['n_all'] = int(len(sub))
    r['share_last_third'] = float((sub.third == 2).mean()); r['rain16_45_median'] = float(sub.rain16_45.median())
    r['events_by_region'] = sub.ar6.value_counts().to_dict(); r['m3_draws'] = min(a.boot, kw.get('m3boot', a.boot)); res.append(r)
    print(tag, '| events', r['n_events'], '| T50', np.round(r['T50'], 2).tolist(), '| dT50', round(r['dT50'], 2), np.round(r.get('dT50_ci', [np.nan] * 2), 2).tolist(),
          '| rate ratio', round(r['M3_rr_wet'], 3), np.round(r.get('M3_ci', [np.nan] * 2), 3).tolist(), '| level', round(r['level_8_15'], 3), flush=True)
    json.dump(res, open(f'{a.out}/{a.name or a.group}.json', 'w'), indent=1)
