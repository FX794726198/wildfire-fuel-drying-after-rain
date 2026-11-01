# Candidate stations from isd-lite file sizes: >= 12 kB (roughly >= 4 obs/day) in >= 38 of 44 years, >= 18 early and >= 19 late years.
import numpy as np, pandas as pd, glob
D = '.'
L = []
for f in sorted(glob.glob(f'{D}/isd/lists/*.tsv')):
    y = int(f.split('/')[-1][:4]); t = pd.read_csv(f, sep='\t', names=['k', 'size'], dtype={'k': str}); t['y'] = y; L.append(t)
L = pd.concat(L); L['sid'] = L.k.str.rsplit('-', n=1).str[0]
g = L[L['size'] >= int(__import__('os').environ.get('MINSIZE', 12000))]
S = g.groupby('sid')['y'].agg(ny='count', n_e=lambda y: (y < 2001).sum(), n_l=lambda y: (y >= 2001).sum()).reset_index()
S = S[(S['ny'] >= 38) & (S['n_e'] >= 18) & (S['n_l'] >= 19)]
S['USAF'], S['WBAN'] = S.sid.str[:6], S.sid.str[7:]
H = pd.read_csv(f'{D}/isd/isd-history.csv', dtype={'USAF': str, 'WBAN': str})
H = H.rename(columns={'LAT': 'lat', 'LON': 'lon', 'ELEV(M)': 'elev', 'CTRY': 'ctry', 'STATION NAME': 'name'})
S = S.merge(H[['USAF', 'WBAN', 'name', 'ctry', 'lat', 'lon', 'elev']], on=['USAF', 'WBAN'], how='left').dropna(subset=['lat', 'lon'])
S = S[(S.lat.abs() > 0.01) | (S.lon.abs() > 0.01)]
C = np.load(f'{D}/era5g/cells.npz'); key = {(int(a), int(o)): k for k, (a, o) in enumerate(zip(C['ia'], C['io']))}
ia = np.rint((S.lat + 90) / 1.5).astype(int); io = (np.rint((S.lon % 360) / 1.5).astype(int)) % 240
S['cell'] = [key.get((a, o), -1) for a, o in zip(ia, io)]
fire = pd.read_csv(f'{D}/out_era5/cells_obs.csv').cell.values
S = S[S.cell.isin(fire)].sort_values('ny', ascending=False).drop_duplicates(subset=['lat', 'lon'])
S['clat'] = C['lat'][S.cell]; S['clon'] = C['lon'][S.cell]
S.to_csv(f'{D}/isd/stations_cand.csv', index=False)
S['zone'] = pd.cut(S.lat, [-60, -23.5, 23.5, 90], labels=['S extratropics', 'tropics', 'N extratropics'])
print('stations', len(S), 'cells', S.cell.nunique(), 'of', len(fire)); print(S.zone.value_counts().to_dict())
print(S.ctry.value_counts().head(40).to_dict())
