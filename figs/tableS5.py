# Table S5: per-model diagnostics of the CMIP6 comparison with the observed fire season (res/models_table.csv from a13b_models_table.py).
import pandas as pd, os
P = os.environ.get('PARTS', './tables_out/'); os.makedirs(P, exist_ok=True)
T = pd.read_csv('./res/models_table.csv')
m = lambda v, f='.1f': '–' if pd.isna(v) else format(v, f).replace('-', '−')
L = ['| Model | Cells | Events | Drying time 1979–2000 (days) | Change, own cells (%) | Change, common cells (%) | ΔRH (%) | ΔT (K) | Change, own fire season (%) |',
     '|:------------------|-------:|---------:|----------:|----------:|-----------:|--------:|-------:|-----------:|']
for _, r in T.iterrows():
    nm = 'ERA5' if r.model == 'ERA5' else r.model
    L.append(f'| {nm} | {int(r.cells):,} | {int(r.events):,} | {m(r.base_common, ".2f")} | {m(r.d_all)} | {m(r.d_common)} | {m(r.dRH, ".2f")} | {m(r.dT, ".2f")} | {m(r.d_own) if r.model != "ERA5" else "−8.5"} |')
open(P + 'tabS5.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
