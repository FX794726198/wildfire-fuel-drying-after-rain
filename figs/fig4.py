# Fig. 4 Human influence, climate models and projections of the drying time after rain (v1.2: DAMIP on fixed events with corrected
# single-forcing responses; new panel d compares stations, ERA5 and CMIP6 on the same cells, season and rules; 3 degC map moved to Fig. S8).
import numpy as np, pandas as pd, sys, os, json, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import plt, MM, OI, HALO, binned, overlap_check
from figmaps import ROB, ASPECT, cellmap, hcolorbar, stipple
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/data/'
m = lambda v, f='.1f': format(v, f).replace('-', '−')
S = json.load(open(D + 'fig4_summary.json')); H = pd.read_csv(D + 'fig4_hist.csv'); H = H[H.use]
G = pd.read_csv(D + 'fig4_gwl_models.csv', index_col=0); F5 = pd.read_csv(D + 'fig4_damip.csv')
SC = json.load(open(D + 'fig4_samecov.json')); SM = pd.read_csv(D + 'fig4_samecov_models.csv', index_col=0)
DEC = pd.read_csv(D + 'fig3_decomp.csv'); neither = DEC[(DEC.rain == 'ERA5 rain') & (DEC.case == 'Neither')].iloc[0]
Pi = pd.read_csv(D + 'fig4_pi.csv') if os.path.exists(D + 'fig4_pi.csv') else None
LE = pd.read_csv(D + 'fig4_le_members.csv')
LEC = {'CanESM5': OI['sky'], 'MIROC6': OI['green'], 'IPSL-CM6A-LR': OI['purple'], 'MPI-ESM1-2-LR': OI['orange']}
WIN = os.environ.get('WIN', 'gwl3.0')
WMM, HMM = 183, 124
fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
top = 73; ht = 47
a = ax_(3, top, 50, ht); b = ax_(72, top, 46, ht); c = ax_(128, top, 54, ht)
d = ax_(50, 9, 60, 53); e = ax_(129, 9, 52, 49)
# ---- a: DAMIP (same rain events; responses of historical - hist-nat, hist-GHG and hist-aer relative to 1850-1900)
rows = [('ERA5 weather', 'obs', OI['verm']), ('Human T and RH response removed', 'ant', 'k'), ('Human warming removed', 'ant_T', '#6E6E6E'),
        ('Human humidity change removed', 'ant_H', '#9E9E9E'), ('Greenhouse-gas response removed', 'ghg', OI['purple']),
        ('Aerosol response removed', 'aer', OI['blue'])]
x0 = -12.5
for i, (lab, kind, cc) in enumerate(rows):
    y = len(rows) - 1 - i
    if kind == 'obs': v = F5.med_ref.iloc[0]; pts = []
    else: k = F5[F5.kind == kind]; v = k.med.median(); pts = k.med.values
    a.barh(y, v, height=0.45, color=cc, alpha=0.85)
    if len(pts): a.scatter(pts, y + np.linspace(-0.13, 0.13, len(pts)), s=5, color='white', edgecolor='k', linewidth=0.4, zorder=3)
    a.text(x0 + 0.2, y + 0.3, lab, fontsize=6.5, ha='left', va='bottom', color='k')
    a.text(0.35, y, m(v), fontsize=6.5, ha='left', va='center', color='k')
a.axvline(0, color='k', lw=0.5); a.set_xlim(x0, 2.6); a.set_ylim(-0.5, len(rows) - 0.1)
a.set_yticks([]); a.spines['left'].set_visible(False)
a.set_xticks([-10, -5, 0]); a.set_xticklabels(['−10', '−5', '0']); a.set_xlabel('Median change in drying time (%)', labelpad=1)
# ---- b: models vs observations, 1979-2000 to 2001-2022
b.axhline(0, color='k', lw=0.4); b.axvline(0, color='k', lw=0.4)
sl, ic = np.polyfit(H.dRH, H.d_med, 1); xx = np.array([S['obs_dRH'], H.dRH.max()])
b.plot(xx, ic + sl * xx, color='#6E6E6E', lw=0.7, ls='--', zorder=2)
for mdl, g in LE.groupby('model'):
    b.scatter(g.dRH, g.d_med, s=4, color=LEC.get(mdl, '#BDBDBD'), alpha=0.75, linewidths=0, zorder=2.5)
b.scatter(H.dRH, H.d_med, s=14, color='#9E9E9E', edgecolor='k', linewidth=0.3, zorder=3)
b.scatter([S['obs_dRH']], [S['obs_d']], s=60, marker='*', color=OI['verm'], edgecolor='k', linewidth=0.3, zorder=4)
b.text(S['obs_dRH'], S['obs_d'] + 0.7, 'ERA5', color=OI['verm'], fontsize=6.5, va='bottom', ha='center')
b.text(0.97, 0.04, f"{S['n_models']} CMIP6 models (grey),\nown fire seasons, $r$ = {S['r_d_dRH']:.2f}", color='#6E6E6E', fontsize=6.5, transform=b.transAxes, ha='right', va='bottom', path_effects=HALO)
b.text(0.03, 0.97, 'Large-ensemble members:', color='#4D4D4D', fontsize=6, transform=b.transAxes, ha='left', va='top')
for i, (mdl, cc) in enumerate(LEC.items()):
    nm = (LE.model == mdl).sum()
    if nm: b.text(0.03, 0.97 - (i + 1) * 0.06, f'{mdl} ({nm})', color=cc, fontsize=6, transform=b.transAxes, ha='left', va='top')
b.set_xlim(-2.6, max(1.0, LE.dRH.max() + 0.2)); b.set_ylim(-12, 6); b.set_xticks([-2, -1, 0, 1]); b.set_xticklabels(['−2', '−1', '0', '1'])
b.set_yticks([-10, -5, 0]); b.set_yticklabels(['−10', '−5', '0'])
b.set_xlabel('Change in daily min. RH (%)', labelpad=1); b.set_ylabel('Change in drying time (%)', labelpad=1)
# ---- c: internal variability in control simulations
if Pi is not None:
    v = Pi.d_med.values; bins = np.arange(-12, 8.1, 1)
    c.hist(np.clip(v, -11.9, 7.9), bins=bins, color='#BDBDBD', edgecolor='white', linewidth=0.3, density=True)
    c.axvline(S['obs_d'], color=OI['verm'], lw=1.0); c.axvline(S['hist_d_med'], color='k', lw=0.8, ls='--')
    c.text(S['obs_d'] + 0.3, 0.97, 'ERA5', color=OI['verm'], fontsize=6.5, transform=c.get_xaxis_transform(), ha='left', va='top')
    c.text(S['hist_d_med'] - 0.3, 0.80, 'Models,\nhistorical', color='k', fontsize=6.5, transform=c.get_xaxis_transform(), ha='right', va='top')
    c.text(0.99, 0.97, f"{S['pi_models']} models,\n{S['pi_windows']:,}\ncontrol windows", fontsize=6.5, color='#6E6E6E', transform=c.transAxes, ha='right', va='top')
c.set_xlim(-12, 8); c.set_xticks([-10, -5, 0, 5]); c.set_xticklabels(['−10', '−5', '0', '5']); c.set_yticks([])
c.spines['left'].set_visible(False); c.set_xlabel('Change between 22-year halves (%)', labelpad=1)
# ---- d: same cells, observed fire-season months and rules for stations, ERA5 and 21 CMIP6 models
GR = [('common', f"All common cells ({SC['common']['cells']:,})", [('era5', 'ERA5', 'era'), ('era5_hsc_0.5', 'ERA5, humidity trend halved', 'eraL'), ('era5_hsc_0', 'ERA5, no humidity trend', 'eraL'), ('models', 'Climate models', 'mod')]),
      ('stationsA', f"Station cells, gaps filled with ERA5 ({SC['stationsA']['cells']})", [('station', 'Stations', 'st'), ('era5_sthours', 'ERA5, station hours', 'eraO'), ('era5', 'ERA5, cells', 'era'), ('models', 'Climate models', 'mod')]),
      ('stationsB', f"Station cells, station only ({SC['stationsB']['cells']})", [('station', 'Stations', 'st'), ('era5_sthours', 'ERA5, same events', 'eraO'), ('models', 'Climate models', 'mod')])]
yy = 0; ylab = []
for g, title, items in GR[::-1]:
    for key, lab, sty in items[::-1]:
        q = SC[g]
        if sty == 'mod':
            v = SM[g].dropna().values; d.scatter(v, yy + np.zeros(len(v)), s=6, color='#BDBDBD', edgecolor='#6E6E6E', linewidth=0.3, zorder=2)
            d.plot([np.median(v)] * 2, [yy - 0.3, yy + 0.3], color='k', lw=1.2, zorder=3)
        else:
            col = {'era': OI['verm'], 'eraL': OI['verm'], 'eraO': '#595959', 'st': OI['purple']}[sty]
            fc = 'white' if sty in ('eraL', 'eraO') else col; ci = q[key + '_ci']
            d.plot(ci, [yy, yy], color=col, lw=0.8, solid_capstyle='butt', alpha=0.6 if sty == 'eraL' else 1)
            d.scatter(q[key], yy, s=14, color=fc, edgecolor=col, linewidth=0.7, zorder=3)
        d.text(-0.02, yy, lab, fontsize=6.5, ha='right', va='center', transform=d.get_yaxis_transform()); yy += 1
    d.text(-0.02, yy - 0.05, title, fontsize=6.5, fontweight='bold', ha='right', va='center', transform=d.get_yaxis_transform()); yy += 1.25
d.axvline(0, color='k', lw=0.5); d.set_ylim(-0.6, yy - 0.6); d.set_yticks([]); d.spines['left'].set_visible(False)
d.set_xlim(-14, 4); d.set_xticks([-10, -5, 0]); d.set_xticklabels(['−10', '−5', '0'])
d.set_xlabel('Median change in drying time, same cells and rules (%)', labelpad=1)
# ---- e: projections by global warming level
lv = [1.0, 1.5, 2.0, 3.0, 4.0]
for j, L in enumerate(lv):
    col = f'gwl_{L:.1f}'; vv = G[col].dropna().values if col in G else np.array([])
    if len(vv) == 0: continue
    q = np.percentile(vv, [10, 25, 50, 75, 90])
    e.plot([L, L], [q[0], q[4]], color='#6E6E6E', lw=0.7); e.add_patch(matplotlib.patches.Rectangle((L - 0.12, q[1]), 0.24, q[3] - q[1], color=OI['blue'], alpha=0.35, lw=0))
    e.plot([L - 0.12, L + 0.12], [q[2], q[2]], color=OI['blue'], lw=1.2)
    e.text(L, 1.2, f'{len(vv)}', fontsize=6.5, color='#6E6E6E', ha='center', va='bottom')
# observation-anchored: observed change plus each model's further change beyond 1 degC
ST = json.load(open(D + 'figS5_txrn.json'))['B_all']['st'][0]                      # station-only change (daily max T / min RH)
for L in lv[1:]:
    col = f'gwl_{L:.1f}'; ok = G[[col, 'gwl_1.0']].dropna()
    for ref, cc, dx in [(S['obs_d'], OI['verm'], 0.13), (ST, OI['purple'], 0.25)]:
        an = 100 * ((1 + ref / 100) * (1 + ok[col] / 100) / (1 + ok['gwl_1.0'] / 100) - 1); q = np.percentile(an, [10, 50, 90])
        e.plot([L + dx, L + dx], [q[0], q[2]], color=cc, lw=0.7); e.scatter([L + dx], [q[1]], s=14, marker='D', color='white', edgecolor=cc, linewidth=0.8, zorder=4)
OCI = json.load(open(D + 'fig4_obs_ci.json'))['obs_d_ci']; STCI = json.load(open(D + 'figS5_txrn.json'))['B_all']['st'][1:]   # 95% intervals
e.plot([0.82, 0.82], OCI, color=OI['verm'], lw=0.8); e.scatter([0.82], [S['obs_d']], s=60, marker='*', color=OI['verm'], edgecolor='k', linewidth=0.3, zorder=4)
e.text(0.97, S['obs_d'], 'ERA5', color=OI['verm'], fontsize=6.5, ha='left', va='center')
e.text(0.04, 0.04, 'Stars: 2001–2022\nBoxes: models\nDiamonds: ERA5 or station change\n+ model change beyond 1 °C', color='#4D4D4D', fontsize=6, transform=e.transAxes, ha='left', va='bottom')
e.plot([0.56, 0.56], STCI, color=OI['purple'], lw=0.8); e.scatter([0.56], [ST], s=30, marker='*', color=OI['purple'], edgecolor='k', linewidth=0.3, zorder=4)
e.text(0.33, ST, 'Stations', color=OI['purple'], fontsize=6.5, ha='center', va='center', rotation=90)
e.axhline(0, color='k', lw=0.4); e.set_xlim(0.2, 4.5); e.set_xticks(lv); e.set_xticklabels(['1', '1.5', '2', '3', '4'])
e.set_ylim(-27, 4); e.set_yticks([-20, -10, 0]); e.set_yticklabels(['−20', '−10', '0'])
e.set_xlabel('Global warming above 1850–1900 (°C)', labelpad=1); e.set_ylabel('Change in drying time\nfrom 1979–2000 (%)', labelpad=1)
for s, xx, yy in [('a', 0.5, HMM - 0.5), ('b', 59, HMM - 0.5), ('c', 120, HMM - 0.5), ('d', 0.5, 66), ('e', 116, 66)]:
    fig.text(xx / WMM, yy / HMM, s, fontsize=9, fontweight='bold', va='top', ha='left')
bad = overlap_check(fig); print('overlap check:', bad or 'none')
fig.savefig(HERE + '/Fig4.png', dpi=600); fig.savefig(HERE + '/Fig4.pdf')
from PIL import Image; im = Image.open(HERE + '/Fig4.png'); im.resize((im.size[0] // 3, im.size[1] // 3)).save(HERE + '/Fig4_prev.png')
