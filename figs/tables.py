# Supplementary tables S1 (alternative inputs and subsets), S2 (trend removal and pre-rain dryness) and S3 (data sets) as Markdown pipe
# tables for the NC manuscript; Table S4 (climate models) is written by tableS4.py once the CMIP6 results are in.
import numpy as np, pandas as pd, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import block_ci
B = './'; O = B + 'out_era5/'; P = os.environ.get('PARTS', './tables_out/'); os.makedirs(P, exist_ok=True)
C = pd.read_csv(B + 'cells_all.csv'); SEL = set(pd.read_csv(B + 'cells_sel.csv').cell)
m = lambda v, f='.1f': format(v, f).replace('-', '−')
def row(name, R, col='rel', note=''):
    lo, hi = block_ci(R, col=col, nb=1000)
    return f'| {name} | {len(R):,} | {m(R[col].median())} ({m(lo)} to {m(hi)}) | {(R[col] < 0).mean() * 100:.0f} |{note}'
def cells(tag, sel=True):
    R = pd.read_csv(O + f'cells_{tag}.csv').merge(C[['cell', 'lat', 'lon']], on='cell')
    if sel: R = R[R.cell.isin(SEL)]
    R['rel'] = 100 * R.d; return R
def from_events(E, col='rec', ref=None):
    E = E.assign(per=np.where(E.y < 2001, 'e', 'l'))
    k = E.groupby(['cell', 'per'])[col].agg(['mean', 'size']).unstack()
    k = k[(k[('size', 'e')] >= 10) & (k[('size', 'l')] >= 10)]
    R = pd.DataFrame(dict(cell=k.index, rel=100 * (k[('mean', 'l')] / k[('mean', 'e')] - 1).values)).merge(C[['cell', 'lat', 'lon']], on='cell'); return R
L = ['| Version | Cells or stations | Median change in drying time, % (95% CI) | Shorter (%) |', '|:------------------------------|------:|:-------------------|------:|']
ref = cells('obs'); L.append(row('Reference: ERA5 weather and rain, 1979–2000 vs 2001–2022', ref))
ch_all = cells('rain_chirps', sel=False); L.append(row('CHIRPS rain, all cells that qualify (1981–2000 vs 2001–2022)', ch_all))
ch = cells('rain_chirps'); L.append(row('CHIRPS rain, cells of the reference set', ch))
L.append(row('ERA5 rain, the same cells', ref[ref.cell.isin(ch.cell)]))
if os.path.exists(O + 'cells_rain_cpc.csv'):
    cp = cells('rain_cpc'); L.append(row('CPC gauge rain, cells of the reference set', cp)); L.append(row('ERA5 rain, the same cells', ref[ref.cell.isin(cp.cell)]))
L.append(row('Daily ERA5 noon wind', cells('windobs')))
L.append(row('Daily maximum temperature and minimum relative humidity', cells('txrn')))
L.append(row('Dryness threshold FFMC 90', cells('rt90')))
E = pd.read_csv(O + 'events_obs.csv.gz'); E = E[E.cell.isin(SEL)]
L.append(row('Events without a further rain ≥ 10 mm within 10 days (ERA5 rain)', from_events(E[~E.further10])))
EC = pd.read_csv(O + 'events_rain_chirps.csv.gz'); EC = EC[EC.cell.isin(SEL)]
L.append(row('Events without a further rain ≥ 10 mm within 10 days (CHIRPS rain)', from_events(EC[~EC.further10])))
L.append(row('FFMC run forward from the rain day without any further rain', from_events(E, 'rec_nr')))
NE = pd.read_csv(B + 'res/nelson_events.csv.gz')
for v, t, nm in [('m1_14', 0.08, '1-h sticks, 8%'), ('m1_14', 0.10, '1-h sticks, 10%'), ('m1_14', 0.12, '1-h sticks, 12%'),
                 ('m10_14', 0.10, '10-h sticks, 10%'), ('m10_14', 0.12, '10-h sticks, 12%'), ('m10_14', 0.15, '10-h sticks, 15%')]:
    L.append(row(f'Nelson model, {nm} moisture at 14:00', from_events(NE[(NE['var'] == v) & (NE.thr == t)], 'rc')))
ST = pd.read_csv('./res/robustness_versions.csv').set_index('version')
for k, nm in [('ERA5', 'Weather stations in seven extratropical regions, station rain and ERA5 humidity (1979–2001 vs 2002–2024)'), ('Station humidity (ISD)', 'The same regions, stations with their own rain and humidity measurements')]:
    r = ST.loc[k]; L.append(f'| {nm} | {int(r.n):,} | {m(r.med)} ({m(r.lo)} to {m(r.hi)}) | {r.shorter * 100:.0f} |')
import json
J = json.load(open(os.path.dirname(os.path.abspath(__file__)) + '/data/figS5_summary.json'))
for t, rn in [('main', 'ERA5 rain'), ('cpc', 'CPC gauge rain')]:
    a = J[t]['all']
    for key, nm in [('st', f'Global ISD stations, measured noon temperature and humidity, {rn}'), ('era', f'ERA5 temperature and humidity at the same stations and hours, {rn}')]:
        v = a[key]; L.append(f"| {nm} | {a['n']:,} | {m(v[0])} ({m(v[1])} to {m(v[2])}) | {a['shorter_' + key] * 100:.0f} |")
open(P + 'tabS1.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
# ---- Table S2: trend removal and pre-rain dryness
L2 = ['| Analysis | Version | Cells or events | Change in drying time, % (95% CI) |', '|:------------------|:---------------------------|------:|:-------------------|']
for rain, tags in [('ERA5 rain', ['obs', 'warm', 'dry', 'neither']), ('CHIRPS rain', ['rain_chirps', 'rain_chirps_warm', 'rain_chirps_dry', 'rain_chirps_neither'])]:
    for lab, t in zip(['Observed weather', 'Humidity trend removed', 'Temperature trend removed', 'Both trends removed'], tags):
        R = cells(t); lo, hi = block_ci(R, nb=1000)
        L2.append(f'| Trend removal, {rain} | {lab} | {len(R):,} cells | {m(R.rel.median())} ({m(lo)} to {m(hi)}) |')
E['fc'] = pd.cut(E.fpre, [85, 88, 91, 94, 200], right=False, labels=['85–88', '88–91', '91–94', '≥ 94'])
rng = np.random.default_rng(5); cl = E.cell.unique(); G = {c: g for c, g in E.groupby('cell')}
def pooled(e):
    k = e.groupby(['fc', e.y >= 2001], observed=True).rec.mean().unstack(); return 100 * (k[True] / k[False] - 1)
base = pooled(E); bs = np.array([pooled(pd.concat([G[c] for c in rng.choice(cl, len(cl))])).reindex(base.index).values for _ in range(200)])
for j, (fc, v) in enumerate(base.items()):
    n = (E.fc == fc).sum(); lo, hi = np.percentile(bs[:, j], [2.5, 97.5])
    L2.append(f'| Pre-rain FFMC class (pooled events) | FFMC {fc} on the day before the rain | {n:,} events | {m(v)} ({m(lo)} to {m(hi)}) |')
open(P + 'tabS2.md', 'w').write('\n'.join(L2) + '\n'); print('\n'.join(L2))
