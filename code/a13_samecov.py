# Same-coverage comparison of the change in drying time (1979-2000 -> 2001-2022) between weather stations, ERA5 and CMIP6, with the
# same inputs (daily maximum temperature and minimum relative humidity, fixed 12 km/h wind), the same observed fire-season months,
# the same cells, drying times capped at 15 days and pooled over each period's events, and the median over cells.
# Cell sets: 'common' = cells where ERA5 and every model have >= 10 events in each period; 'stations' = common cells with a valid
# station (version A, ERA5-assisted; version B, station-only). Output res/samecov_*.csv/json.
import numpy as np, pandas as pd, glob, os, sys, json
sys.path.insert(0, './figs'); from figcommon import block_ci
D = '.'; CS = pd.read_csv(f'{D}/code/cells_sel.csv'); LL = CS.set_index('cell')[['lat', 'lon']]
def pooled(cell, y, s, n):
    t = pd.DataFrame(dict(cell=cell, late=y >= 2001, s=s, n=n)).groupby(['cell', 'late'])[['s', 'n']].sum().unstack()
    t = t[(t[('n', False)] >= 10) & (t[('n', True)] >= 10)]
    return (t[('s', True)] / t[('n', True)]) / (t[('s', False)] / t[('n', False)]) * 100 - 100
# ERA5 hourly Tmax/RHmin with the observed fire season (fixed-event runs: obs and humidity-trend scalings)
ev = pd.read_csv(f'{D}/fx/events_txrn.csv.gz')
E = {}
for v in ['obs', 'hsc_0.5', 'hsc_0', 'hsc_0.5_far', 'hsc_0_far']:
    r = np.minimum(np.load(f'{D}/fx/rec_txrn_{v}.npy').astype(float), 15); E[v] = pooled(ev.cell.values, ev.y.values, r, np.ones(len(r)))
# CMIP6 models with the observed fire season
M = {}; MH = {}; MT = {}
for f in sorted(glob.glob(f'{D}/fs/*_hist*_*.csv.gz')):
    mdl = os.path.basename(f).split('_')[0]; g0 = pd.read_csv(f); g0 = g0[(g0.y >= 1979) & (g0.y <= 2022)]; g = g0[g0.n.notna()]
    M[mdl] = pooled(CS.cell.values[g.cell.values], g.y.values, (g.rec15 * g.n).values, g.n.values)
    h = g0.assign(cell=CS.cell.values[g0.cell.values], late=g0.y >= 2001).groupby(['cell', 'late'])[['H', 'T']].mean().unstack()
    MH[mdl] = h[('H', True)] - h[('H', False)]; MT[mdl] = h[('T', True)] - h[('T', False)]   # fire-season daily min RH (pp), max T (K)
common = set(E['obs'].index)
for v in M.values(): common &= set(v.index)
common = sorted(common)
# stations
R = pd.read_csv(f'{D}/res/isd_txrn.csv')
okA = R.dA_st.notna() & R.dA_era.notna() & (R.cov_e >= 0.7) & (R.cov_l >= 0.7); okB = R.dB_st.notna() & R.dB_era.notna()
SA = R[okA].groupby('cell')[['dA_st', 'dA_era']].mean() * 100; SB = R[okB].groupby('cell')[['dB_st', 'dB_era']].mean() * 100
out = {}; rows = []
for name, cells in [('common', common), ('stationsA', sorted(set(common) & set(SA.index))), ('stationsB', sorted(set(common) & set(SB.index)))]:
    X = pd.DataFrame(index=cells).join(LL); X['era5'] = E['obs'].reindex(cells)
    for v in ['hsc_0.5', 'hsc_0', 'hsc_0.5_far', 'hsc_0_far']: X['era5_' + v] = E[v].reindex(cells)
    if name == 'stationsA': X['station'] = SA.dA_st.reindex(cells); X['era5_sthours'] = SA.dA_era.reindex(cells)
    if name == 'stationsB': X['station'] = SB.dB_st.reindex(cells); X['era5_sthours'] = SB.dB_era.reindex(cells)
    for m, v in M.items(): X['m_' + m] = v.reindex(cells)
    mods = [c for c in X if c.startswith('m_')]; med = X[mods].median()
    res = dict(cells=len(cells), models=len(mods), model_med=float(med.median()), model_min=float(med.min()), model_max=float(med.max()))
    res['model_med_noKIOST'] = float(med.drop('m_KIOST-ESM', errors='ignore').median())   # KIOST-ESM's second period ends in 2017
    for c in [c for c in X if not c.startswith('m_') and c not in ('lat', 'lon')]:
        res[c] = float(X[c].median()); res[c + '_ci'] = block_ci(X.dropna(subset=[c]), c, nb=1000).tolist()
        res[c + '_beyond_models'] = int((med < X[c].median()).sum())          # models with a more negative median
    res['model_dRH'] = float(pd.DataFrame({m: v.reindex(cells) for m, v in MH.items()}).median().median())
    res['model_dT'] = float(pd.DataFrame({m: v.reindex(cells) for m, v in MT.items()}).median().median())
    if name == 'stationsB':
        # T/RH changes from the station-only values on valid fire-season days, ERA5 on the same days (v1.3)
        for a, b in [('station_dRH', 'dRH_stB'), ('era5_dRH', 'dRH_eraB'), ('station_dT', 'dT_stB'), ('era5_dT', 'dT_eraB')]: res[a] = float(R[okB].groupby('cell')[b].mean().reindex(cells).median())
    if 'station' in X:
        X['pd'] = X.station - X.era5_sthours; res['paired_st_era'] = float(X.pd.median()); res['paired_ci'] = block_ci(X, 'pd', nb=1000).tolist()
    out[name] = res; X.to_csv(f'{D}/res/samecov_{name}.csv'); rows.append(med.rename(name))
pd.concat(rows, axis=1).to_csv(f'{D}/res/samecov_models.csv'); json.dump(out, open(f'{D}/res/samecov.json', 'w'), indent=1)
for k, v in out.items(): print(k, {a: (round(b, 2) if isinstance(b, float) else b) for a, b in v.items() if not a.endswith('_ci')})
