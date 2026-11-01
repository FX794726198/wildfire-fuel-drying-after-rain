# Supplementary figures S7 (CMIP6 models against observations, 1979-2022) and S8 (projected change at four warming levels).
import numpy as np, pandas as pd, sys, os, json, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import plt, MM, OI, HALO, binned, overlap_check
from figmaps import ROB, ASPECT, cellmap, hcolorbar, stipple
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/data/'
m = lambda v, f='.1f': format(v, f).replace('-', '−')
S = json.load(open(D + 'fig4_summary.json')); H = pd.read_csv(D + 'fig4_hist.csv'); H = H[H.use].sort_values('d_med').reset_index(drop=True)
SE = pd.read_csv(D + 'figS6_series.csv', index_col=0); MP = pd.read_csv(D + 'fig4_map.csv')
def save(fig, name):
    bad = overlap_check(fig); print(name, 'overlap check:', bad or 'none')
    fig.savefig(f'{HERE}/Fig{name}.png', dpi=600); fig.savefig(f'{HERE}/Fig{name}.pdf')
    from PIL import Image; im = Image.open(f'{HERE}/Fig{name}.png'); im.resize((im.size[0] // 3, im.size[1] // 3)).save(f'{HERE}/Fig{name}_prev.png')
# ---------------- S6
WMM, HMM = 183, 66
fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
a = ax_(13, 21, 62, 42); b = ax_(92, 11, 38, 52); c = ax_(144, 11, 38, 52)
y = np.arange(len(H))
a.bar(y, H.d_med, color='#9E9E9E', width=0.7); a.scatter(y, H.d_med, s=5, color='k', zorder=3, linewidths=0)
a.axhline(S['obs_d'], color=OI['verm'], lw=0.9); a.axhline(S['hist_d_med'], color='k', lw=0.7, ls='--'); a.axhline(0, color='k', lw=0.4)
a.set_xticks(y); a.set_xticklabels(H.model, rotation=90, fontsize=6); a.set_xlim(-0.7, len(H) - 0.3); a.tick_params(axis='x', pad=1)
a.set_ylim(-12, 4); a.set_yticks([-10, -5, 0]); a.set_yticklabels(['−10', '−5', '0']); a.set_ylabel('Change in drying time,\n1979–2000 to 2001–2022 (%)', labelpad=1)
a.text(0.98, 0.08, 'Observed (ERA5)', color=OI['verm'], transform=a.transAxes, fontsize=6.5, ha='right', va='bottom')
a.text(0.98, S['hist_d_med'] - 0.4, f"Model median {m(S['hist_d_med'])}%", color='k', transform=a.get_yaxis_transform(), fontsize=6.5, ha='right', va='top')
for ax, col, ov, lb, lim, rr in [(b, 'dT', S['obs_dT'], 'Change in daily max. temperature (K)', (0, 1.6), S['r_d_dT']),
                                 (c, 'dRH', S['obs_dRH'], 'Change in daily min. RH (%)', (-2.5, 1.5), S['r_d_dRH'])]:
    ax.scatter(H[col], H.d_med, s=12, color='#9E9E9E', edgecolor='k', linewidth=0.3, zorder=3)
    ax.scatter([ov], [S['obs_d']], s=50, marker='*', color=OI['verm'], edgecolor='k', linewidth=0.3, zorder=4)
    ax.text(ov, S['obs_d'] + 0.7, 'ERA5', color=OI['verm'], fontsize=6.5, ha='center', va='bottom')
    ax.text(0.97, 0.04, f'$r$ = {m(rr, ".2f")}', color='#6E6E6E', fontsize=6.5, transform=ax.transAxes, ha='right', va='bottom')
    ax.set_xlim(*lim); ax.set_ylim(-12, 4); ax.set_yticks([-10, -5, 0]); ax.set_yticklabels(['−10', '−5', '0'])
    ax.set_xlabel(lb, labelpad=1); ax.axhline(0, color='k', lw=0.4)
b.set_ylabel('Change in drying time (%)', labelpad=1)
fig.text(0.5 / WMM, (HMM - 0.5) / HMM, 'a', fontsize=9, fontweight='bold', va='top'); fig.text(80 / WMM, (HMM - 0.5) / HMM, 'b', fontsize=9, fontweight='bold', va='top')
fig.text(134 / WMM, (HMM - 0.5) / HMM, 'c', fontsize=9, fontweight='bold', va='top')
save(fig, 'S7')
# ---------------- S7
BOUNDS = [-30, -20, -15, -10, -5, -2, 2, 5, 10, 15, 20, 30]; CM, NORM = binned('BrBG', BOUNDS)
WMM = 183; mw = 90; mh = mw * ASPECT; gap = 6.5; HMM = 2 * (mh + gap) + 12
fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
L = []
for i, w in enumerate(['gwl1.5', 'gwl2.0', 'gwl3.0', 'gwl4.0']):
    r, cidx = divmod(i, 2); x = 1 + cidx * (mw + 1); yy = HMM - (r + 1) * (mh + gap) + 1
    ax = ax_(x, yy, mw, mh, projection=ROB); k = MP[MP.win == w]
    if len(k):
        cellmap(ax, k.lat, k.lon, k.med.values, CM, NORM)
        stipple(ax, k.lat, k.lon, k.neg >= 0.8)
        ax.text(0.06, 0.32, f'Median {m(k.med.median())}%\nDecline: {(k.med < 0).mean() * 100:.0f}%\n≥80% agree: {(k.neg >= 0.8).mean() * 100:.0f}%', transform=ax.transAxes, fontsize=6.5, ha='left', va='top', path_effects=HALO)
    nm = S['gwl'][f'gwl_{w[3:]}'][3]
    ax.text(0.5, 1.02, f'{w[3:].rstrip("0").rstrip(".")} °C above 1850–1900 ({nm} models)', transform=ax.transAxes, ha='center', va='bottom', fontsize=7)
    L.append(('abcd'[i], x - 0.5, min(yy + mh + gap - 0.5, HMM - 0.5)))
hcolorbar(fig, ax_(46, 7.5, 91, 2.2), CM, NORM, BOUNDS, 'Change in drying time from 1979–2000, multi-model median (%)')
for s, x, yy in L: fig.text(x / WMM, yy / HMM, s, fontsize=9, fontweight='bold', va='top', ha='left')
save(fig, 'S8')
