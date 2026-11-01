# CMIP6 daily-FFMC results (from e3_cmip6.py) compared with the observed (ERA5 Tmax + RHmin, 'txrn') wet-fuel-day change.
# 1. historical+SSP: per-model change 1979-2000 -> 2001-2022 (median over the analysis cells), global yearly series and trend,
#    fire-season dRH and dT; model vs observed.
# 2. piControl: distribution of the same statistics in all 44-year windows (step 2 yr) -> detection (signal / noise).
# 3. Projections: change vs 1979-2000 for 2041-2060, 2081-2100 and at 1.5/2/3/4 C global warming levels (Hauser et al. windows).
# Usage: python a5_cmip6.py <dir with *.csv.gz> [tag]
import numpy as np, pandas as pd, glob, os, sys
D = sys.argv[1]; tag = sys.argv[2] if len(sys.argv) > 2 else 'cmip6'
os.makedirs('res', exist_ok=True); SEL = set(pd.read_csv('cells_sel.csv').cell)
def periods(g, a, b):
    """per-cell event-weighted mean rec15 in year ranges a and b, relative change; >=10 events in each."""
    x = g[g.n.notna()].assign(s=lambda q: q.rec15 * q.n)
    def pm(r):
        z = x[(x.y >= r[0]) & (x.y <= r[1])]
        k = z.groupby('cell').agg(s=('s', 'sum'), n=('n', 'sum')); k['r'] = k.s / k.n
        h = g[(g.y >= r[0]) & (g.y <= r[1])].groupby('cell')[['H', 'T']].mean()
        return k.join(h)
    j = pm(a).join(pm(b), lsuffix='_e', rsuffix='_l', how='inner')
    j = j[(j.n_e >= 10) & (j.n_l >= 10)]
    j['d'] = j.r_l / j.r_e - 1; return j
def series(g, ref=(1979, 2022)):
    x = g[g.n.notna()].copy(); m = x[(x.y >= ref[0]) & (x.y <= ref[1])].groupby('cell').rec15.mean()
    x['a'] = x.rec15 / x.cell.map(m) - 1; return x.groupby('y').a.mean()
def trend(s, a=1979, b=2022):
    s = s[(s.index >= a) & (s.index <= b)].dropna(); return np.polyfit(s.index, s.values, 1)[0] * 10 * 100
# ---- observed comparator (ERA5 Tmax/RHmin, same FFMC proxy as the models); skipped when absent (e.g. on DelftBlue)
OBS = None
if os.path.exists('out_era5/events_txrn.csv.gz'):
    E = pd.read_csv('out_era5/events_txrn.csv.gz'); E = E[E.cell.isin(SEL)]
    og = E.groupby(['cell', 'y']).agg(rec15=('rec', lambda x: x.clip(upper=15).mean()), n=('rec', 'size')).reset_index()
    Ro = pd.read_csv('out_era5/cells_txrn.csv'); Ro = Ro[Ro.cell.isin(SEL)]
    og = og.merge(Ro[['cell']], on='cell'); og['H'] = np.nan; og['T'] = np.nan
    po = periods(og, (1979, 2000), (2001, 2022)); so = series(og)
    OBS = dict(model='ERA5', scen='obs', d_med=po.d.median() * 100, shorter=(po.d < 0).mean(), trend=trend(so),
               dRH=(Ro.H_l - Ro.H_e).median(), dT=(Ro.T_l - Ro.T_e).median(), cells=len(po))
    print('OBS', OBS)
rows, ser, pis, proj, cel, pcel = [], {}, [], [], [], []
G = pd.read_csv('gwl/cmip6_warming_levels_all_ens_1850_1900.csv', comment='#', skipinitialspace=True)
for f in sorted(glob.glob(f'{D}/*.csv.gz')):
    mdl, sc, mem = os.path.basename(f)[:-7].split('_')
    g = pd.read_csv(f); g = g[g.cell.isin(SEL)]
    if sc in ('hist245', 'hist585'):
        p = periods(g, (1979, 2000), (2001, 2022)); s = series(g); ser[f'{mdl}_{sc}'] = s
        cel.append(p[['d', 'n_e', 'n_l', 'H_e', 'H_l', 'T_e', 'T_l']].assign(model=mdl, scen=sc).reset_index())
        rows.append(dict(model=mdl, scen=sc, member=mem, d_med=p.d.median() * 100, shorter=(p.d < 0).mean(), trend=trend(s),
                         dRH=(p.H_l - p.H_e).median(), dT=(p.T_l - p.T_e).median(), cells=len(p)))
        out = dict(model=mdl, scen=sc, member=mem)
        for a, b in [(2021, 2040), (2041, 2060), (2081, 2100)]:
            q = periods(g, (1979, 2000), (a, b)); out[f'd_{a}_{b}'] = q.d.median() * 100; out[f'dT_{a}_{b}'] = (q.T_l - q.T_e).median()
            pcel.append(q[['d']].assign(model=mdl, scen=sc, win=f'{a}-{b}').reset_index())
            out[f'dRH_{a}_{b}'] = (q.H_l - q.H_e).median()
        ssp = 'ssp245' if sc == 'hist245' else 'ssp585'
        w = G[(G.model == mdl) & (G.ensemble == mem) & (G.exp == ssp)]
        for _, r in w.iterrows():
            if r.end_year <= g.y.max():
                q = periods(g, (1979, 2000), (int(r.start_year), int(r.end_year))); out[f'gwl_{r.warming_level:.1f}'] = q.d.median() * 100
                pcel.append(q[['d']].assign(model=mdl, scen=sc, win=f'gwl{r.warming_level:.1f}').reset_index())
        proj.append(out); print(mdl, sc, rows[-1], flush=True)
    elif sc == 'picontrol':
        y0, y1 = g.y.min(), g.y.max(); s = series(g, (y0, y1))
        for a in range(y0, y1 - 43, 2):
            p = periods(g, (a, a + 21), (a + 22, a + 43))
            pis.append(dict(model=mdl, start=a, d_med=p.d.median() * 100, trend=trend(s, a, a + 43), cells=len(p)))
        print(mdl, 'piControl windows', len(pis), flush=True)
H = pd.DataFrame(rows + ([OBS] if OBS else []))
H.to_csv(f"res/{tag}_hist.csv", index=False)
pd.DataFrame(ser).to_csv(f'res/{tag}_series.csv'); pd.concat(cel).to_csv(f'res/{tag}_cells.csv.gz', index=False) if cel else None; pd.concat(pcel).to_csv(f'res/{tag}_proj_cells.csv.gz', index=False) if pcel else None; pd.DataFrame(proj).to_csv(f'res/{tag}_proj.csv', index=False)
if pis:
    P = pd.DataFrame(pis); P.to_csv(f'res/{tag}_pi.csv', index=False)
    print('piControl d_med sd', P.d_med.std(), 'trend sd', P.trend.std())
