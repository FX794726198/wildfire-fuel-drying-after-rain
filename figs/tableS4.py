# Table S4: climate models used (DAMIP monthly, CMIP6 daily) and their results, as a Markdown pipe table.
import numpy as np, pandas as pd, os
F = './figs/data/'; P = os.environ.get('PARTS', './tables_out/'); os.makedirs(P, exist_ok=True)
m = lambda v, f='.1f': '–' if v is None or not np.isfinite(v) else format(0.0 if abs(v) < 0.5 * 10 ** -int(f[1]) else v, f).replace('-', '−')
DAMIP = {'CESM2': '10/3/3/–', 'CNRM-CM6-1': '10/10/9/10', 'CanESM5': '10/10/10/10', 'GFDL-ESM4': '3/3/1/1', 'HadGEM3-GC31-LL': '5/4/4/4',
         'IPSL-CM6A-LR': '10/10/10/10', 'MRI-ESM2-0': '10/3/3/3'}
H = pd.read_csv(F + 'fig4_hist.csv'); G = pd.read_csv(F + 'fig4_gwl_models.csv', index_col=0)
Pi = pd.read_csv(F + 'fig4_pi.csv') if os.path.exists(F + 'fig4_pi.csv') else None
models = sorted(set(DAMIP) | set(H.model) | (set(Pi.model) if Pi is not None else set()), key=str.lower)
L = ['| Model | DAMIP members | Daily member | Scenarios | Change (%) | ΔRH (%) | ΔT (K) | 2 °C (%) | 3 °C (%) | Control windows |',
     '|:----------------|:----------|:------------|:----------------------|--------:|--------:|-------:|-------:|-------:|---------:|']
for mdl in models:
    h = H[(H.model == mdl) & H.use]; h = h.iloc[0] if len(h) else None
    sc = ', '.join(s.replace('hist', 'SSP') .replace('SSP245', 'SSP2-4.5').replace('SSP585', 'SSP5-8.5') for s in sorted(H[H.model == mdl].scen.unique()))
    if Pi is not None and mdl in set(Pi.model): sc = (sc + ', ' if sc else '') + 'piControl'
    mem = ', '.join(sorted(H[H.model == mdl].member.unique())) if (H.model == mdl).any() else '–'
    g2 = G.loc[mdl, 'gwl_2.0'] if mdl in G.index and 'gwl_2.0' in G else np.nan; g3 = G.loc[mdl, 'gwl_3.0'] if mdl in G.index and 'gwl_3.0' in G else np.nan
    npi = int((Pi.model == mdl).sum()) if Pi is not None and mdl in set(Pi.model) else None
    L.append(f"| {mdl} | {DAMIP.get(mdl, '–')} | {mem} | {sc or '–'} | {m(h.d_med) if h is not None else '–'} | {m(h.dRH, '.2f') if h is not None else '–'} | "
             f"{m(h.dT, '.2f') if h is not None else '–'} | {m(g2)} | {m(g3)} | {npi if npi else '–'} |")
open(P + 'tabS4.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
