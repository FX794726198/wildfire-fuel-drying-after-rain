# Decompose the change in fire-season dry-fuel days (FFMC >= 85) between 1979-2000 and 2001-2022 into the parts due to
# (a) the number of soaking-rain events, (b) the wet-fuel days per event (faster recovery) and (c) other wet days.
# Usage (from nc/): python code/a7_decomp.py [tag]   (tag = obs by default; reads out_era5/ffmc_<tag>.npz and events_<tag>.csv.gz)
import numpy as np, pandas as pd, sys, json
tag = sys.argv[1] if len(sys.argv) > 1 else 'obs'; O = 'out_era5/'
Z = np.load(O + f'ffmc_{tag}.npz'); F = Z['F']; smon = Z['smon']
d = pd.date_range('1979-01-01', '2022-12-31'); mon = d.month.values; yr = d.year.values
E = pd.read_csv(O + f'events_{tag}.csv.gz'); R = pd.read_csv(O + f'cells_{tag}.csv'); cells = R.cell.values
rows = []
for c in cells:
    s = smon[c, mon - 1]; wet = s & np.isfinite(F[c]) & (F[c] < 85); valid = s & np.isfinite(F[c])
    inev = np.zeros(len(d), bool)
    e = E[E.cell == c]
    for i, r in zip(e.i.values, e.rec.values): inev[i:i + r] = True
    for per, k in [('e', yr <= 2000), ('l', yr >= 2001)]:
        ny = len(np.unique(yr[k]))
        rows.append(dict(cell=c, per=per, season=valid[k].sum() / ny, wet=wet[k].sum() / ny, wet_ev=(wet & inev)[k].sum() / ny,
                         n=((e.y <= 2000) if per == 'e' else (e.y >= 2001)).sum() / ny))
D = pd.DataFrame(rows).pivot(index='cell', columns='per')
out = pd.DataFrame(index=D.index)
for v in ['season', 'wet', 'wet_ev', 'n']: out[v + '_e'] = D[(v, 'e')]; out[v + '_l'] = D[(v, 'l')]
out['oth_e'] = out.wet_e - out.wet_ev_e; out['oth_l'] = out.wet_l - out.wet_ev_l
out['D_e'] = out.wet_ev_e / out.n_e; out['D_l'] = out.wet_ev_l / out.n_l
out['dDry'] = -(out.wet_l - out.wet_e)                                   # change in dry-fuel days per fire season
out['from_n'] = -(out.D_e * (out.n_l - out.n_e))                          # fewer or more soaking rains
out['from_D'] = -(out.n_e * (out.D_l - out.D_e))                          # faster recovery after each rain
out['from_x'] = -((out.n_l - out.n_e) * (out.D_l - out.D_e))              # interaction
out['from_oth'] = -(out.oth_l - out.oth_e)                                # other wet days (light rain, humid spells)
out = out.reset_index(); out.to_csv(f'res/decomp_{tag}.csv', index=False)
w = np.cos(np.deg2rad(pd.read_csv('cells_all.csv').set_index('cell').loc[out.cell, 'lat'].values))
S = {k: float(np.average(out[k], weights=w)) for k in ['season_e', 'wet_e', 'wet_l', 'wet_ev_e', 'wet_ev_l', 'n_e', 'n_l', 'D_e', 'D_l', 'dDry', 'from_n', 'from_D', 'from_x', 'from_oth']}
S['share_D'] = S['from_D'] / S['dDry']; S['share_n'] = S['from_n'] / S['dDry']; S['share_oth'] = S['from_oth'] / S['dDry']
S['median_dDry'] = float(out.dDry.median()); S['median_from_D'] = float(out.from_D.median()); S['frac_cells_from_D_pos'] = float((out.from_D > 0).mean())
json.dump(S, open(f'res/decomp_{tag}.json', 'w'), indent=1); print(json.dumps(S, indent=1))
