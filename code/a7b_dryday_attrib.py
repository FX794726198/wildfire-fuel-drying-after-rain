# Change in fire-season dry-fuel days (FFMC >= 85) between 1979-2000 and 2001-2022 with observed weather and with the
# temperature and/or humidity trends removed (same fire-season months and cells as the observed run).
import numpy as np, pandas as pd, json
O = './out_era5/'
d = pd.date_range('1979-01-01', '2022-12-31'); mon = d.month.values; yr = d.year.values
cells = pd.read_csv(O + 'cells_obs.csv').cell.values
smon = np.load(O + 'ffmc_obs.npz')['smon'][cells]; S = smon[:, mon - 1]
lat = pd.read_csv('./cells_all.csv').set_index('cell').loc[cells, 'lat'].values; w = np.cos(np.deg2rad(lat))
out = {}; per = {}
for v in ['obs', 'neither', 'warm', 'dry']:
    F = np.load(O + f'ffmc_{v}.npz')['F'][cells].astype(np.float32)
    ok = np.isfinite(F) & S; dry = ok & (F >= 85)
    e, l = yr <= 2000, yr >= 2001
    de = dry[:, e].sum(1) / 22; dl = dry[:, l].sum(1) / 22; se = ok[:, e].sum(1) / 22; sl = ok[:, l].sum(1) / 22
    per[v] = dl - de
    out[v] = dict(dry_e=float(np.average(de, weights=w)), dry_l=float(np.average(dl, weights=w)), season=float(np.average(se, weights=w)),
                  dDry=float(np.average(dl - de, weights=w)), median_dDry=float(np.median(dl - de)), frac_up=float((dl > de).mean()))
    del F
for v in ['neither', 'warm', 'dry']:
    out[v]['share_of_obs_removed'] = 1 - out[v]['dDry'] / out['obs']['dDry']
pd.DataFrame({'cell': cells, **{f'dDry_{v}': per[v] for v in per}}).to_csv('./res/dryday_attrib_cells.csv', index=False)
json.dump(out, open('./res/dryday_attrib.json', 'w'), indent=1); print(json.dumps(out, indent=1))
