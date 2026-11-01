# Fig. 3 Warmer, drier air shortens the drying time after rain, and heavier rain cannot restore it.
import numpy as np, pandas as pd, sys, os, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import plt, MM, OI, HALO, binned, overlap_check
import cartopy.crs as ccrs, cartopy.feature as cfeature
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/data/'
R = pd.read_csv(D + 'fig3_cells.csv'); DEC = pd.read_csv(D + 'fig3_decomp.csv'); Q = pd.read_csv(D + 'fig3_scaling.csv')
CU = pd.read_csv(D + 'fig3_rain_curves.csv'); FIT = pd.read_csv(D + 'fig3_rain_fit.csv', index_col=0).iloc[:, 0]
m = lambda v, f='.1f': format(v, f).replace('-', '−')
BOUNDS = [-8, -6, -4, -3, -2, -1, 1, 2, 3, 4]
CM, NORM = binned('BrBG', BOUNDS)
PC = ccrs.PlateCarree(); ROB = ccrs.Robinson(central_longitude=10)
yproj = lambda lat: ROB.transform_points(PC, np.zeros_like(np.asarray(lat, float)), np.asarray(lat, float))[:, 1]
Y0, Y1 = yproj([-57, 80]); X0, X1 = ROB.x_limits
WMM, HMM = 183, 104
fig = plt.figure(figsize=(WMM * MM, HMM * MM))
ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
mw = 117; mh = mw * (Y1 - Y0) / (X1 - X0); my = HMM - 1.5 - mh
ax_map = ax_(1, my, mw, mh, projection=ROB)
ax_b = ax_(127, my + 1, 55, mh - 1)
hb = my - 25
ax_c = ax_(13, 9, 44, hb); ax_d = ax_(72, 9, 52, hb); ax_e = ax_(150, 9, 32, hb)
# ---- a: change in fire-season noon relative humidity
ax_map.set_xlim(X0, X1); ax_map.set_ylim(Y0, Y1); ax_map.set_xticks([]); ax_map.set_yticks([])
ax_map.add_feature(cfeature.LAND.with_scale('110m'), facecolor='#D9D9D9', edgecolor='none', zorder=0)
grid = np.full((121, 240), np.nan)
ia = np.rint((R.lat.values + 90) / 1.5).astype(int); io = np.rint((R.lon.values % 360) / 1.5).astype(int) % 240
grid[ia, io] = R.dRH.values
lon_e = np.arange(241) * 1.5 - 0.75; lat_e = np.clip(np.arange(122) * 1.5 - 90.75, -90, 90)
ax_map.pcolormesh(lon_e, lat_e, grid, cmap=CM, norm=NORM, transform=PC, zorder=1, rasterized=True)
ax_map.coastlines('110m', linewidth=0.3, color='#4D4D4D', zorder=2); ax_map.spines['geo'].set_linewidth(0.5)
ax_map.gridlines(crs=PC, linewidth=0.3, color='#B3B3B3', linestyle='--', xlocs=range(-180, 181, 60), ylocs=[-30, 0, 30, 60], zorder=1.5)
ax_map.text(0.075, 0.36, f'Median {m(R.dRH.median())}%\n{(R.dRH < 0).mean() * 100:.0f}% of cells drier\nWarming {m(R.dT.median(), ".2f")} K',
            transform=ax_map.transAxes, fontsize=6.5, ha='left', va='top', path_effects=HALO)
cax = ax_(22, my - 5.2, 76, 2.2)
cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=NORM, cmap=CM), cax=cax, orientation='horizontal', extend='both', ticks=BOUNDS)
cb.ax.set_xticklabels([m(b, 'd') for b in BOUNDS]); cb.ax.tick_params(labelsize=6.5, length=2, pad=1); cb.outline.set_linewidth(0.4)
cb.set_label('Change in fire-season noon relative humidity, 2001–2022 vs 1979–2000 (%)', fontsize=7, labelpad=1.5)
# ---- b: removing the temperature or humidity trends
cases = ['Observed', 'Warming only', 'Drying only', 'Neither']; lab = {'Observed': 'Observed weather', 'Warming only': 'Humidity trend removed',
                                                                     'Drying only': 'Temperature trend removed', 'Neither': 'Both trends removed'}
ax_b.axvline(0, color='k', lw=0.5)
for j, (rain, c, off) in enumerate([('ERA5 rain', 'k', 0.13), ('CHIRPS rain', OI['blue'], -0.13)]):
    for i, cs in enumerate(cases):
        r = DEC[(DEC.rain == rain) & (DEC.case == cs)].iloc[0]; y = len(cases) - 1 - i + off
        ax_b.plot([r.lo, r.hi], [y, y], color=c, lw=0.9, solid_capstyle='butt'); ax_b.scatter(r.med, y, s=12, color=c, zorder=3, edgecolor='none')
for i, cs in enumerate(cases):
    ax_b.text(-11.8, len(cases) - 1 - i + 0.3, lab[cs], fontsize=6.5, ha='left', va='bottom', path_effects=HALO)
ax_b.set_ylim(-0.45, len(cases) - 0.05); ax_b.set_yticks([]); ax_b.spines['left'].set_visible(False)
ax_b.set_xlim(-12, 4.2); ax_b.set_xticks([-10, -5, 0]); ax_b.set_xticklabels(['−10', '−5', '0'])
ax_b.set_xlabel('Median change in drying\ntime after rain (%)', labelpad=1)
ax_b.text(0.99, 0.995, 'ERA5 rain', color='k', fontsize=6.5, transform=ax_b.transAxes, ha='right', va='top', path_effects=HALO)
ax_b.text(0.99, 0.935, 'CHIRPS rain', color=OI['blue'], fontsize=6.5, transform=ax_b.transAxes, ha='right', va='top', path_effects=HALO)
# ---- c: analytic drying law vs observed (deciles of predicted change)
ax_c.plot([-25, 8], [-25, 8], color='#9E9E9E', lw=0.6, ls='--')
ax_c.errorbar(100 * Q.pred, 100 * Q.obs, yerr=[100 * (Q.obs - Q.lo), 100 * (Q.hi - Q.obs)], fmt='o', ms=3, color=OI['verm'], ecolor=OI['verm'], elinewidth=0.6, capsize=0)
ax_c.set_xlim(-24, 6); ax_c.set_ylim(-24, 10); ax_c.set_xticks([-20, -10, 0]); ax_c.set_yticks([-20, -10, 0, 10])
ax_c.set_xticklabels(['−20', '−10', '0']); ax_c.set_yticklabels(['−20', '−10', '0', '10'])
ax_c.set_xlabel('Predicted by drying relation (%)', labelpad=1); ax_c.set_ylabel('Observed change (%)', labelpad=1)
ax_c.text(0.04, 0.97, '−2.1% per 1% lower RH\n−4.8% per K warmer', transform=ax_c.transAxes, ha='left', va='top', fontsize=6.5, color=OI['verm'])
ax_c.text(0.96, 0.04, 'Deciles of\n1,704 cells', transform=ax_c.transAxes, ha='right', va='bottom', fontsize=6.5, color='#6E6E6E')
# ---- d: drying time against rain amount, early vs late
for late, c, lb in [(0, '#7F7F7F', '1979–2000'), (1, OI['verm'], '2001–2022')]:
    k = CU[(CU.late == late) & (CU.n >= 300)]
    ax_d.fill_between(k.p, k.r - 1.96 * k.se, k.r + 1.96 * k.se, color=c, alpha=0.2, lw=0)
    ax_d.plot(k.p, k.r, color=c, lw=1.0, marker='o', ms=2.5)
    ax_d.text(0.04, 0.985 - 0.085 * late, lb, transform=ax_d.transAxes, color=c, fontsize=6.5, ha='left', va='top')
ax_d.set_xscale('log'); ax_d.set_xticks([10, 20, 50]); ax_d.set_xticklabels(['10', '20', '50']); ax_d.minorticks_off()
ax_d.set_xlim(9.5, 70); ax_d.set_ylim(3.0, 4.6)
ax_d.set_xlabel('Rain amount (mm)', labelpad=1); ax_d.set_ylabel('Drying time after rain (days)', labelpad=1)
e0, l0 = CU[CU.late == 0].iloc[0], CU[(CU.late == 1) & (CU.n >= 300)]
px = np.exp(np.interp(e0.r, l0.r.values, np.log(l0.p.values)))           # rain at which the late curve reaches the early 10-12.5 mm value
ax_d.annotate('', xy=(px, e0.r), xytext=(e0.p, e0.r), arrowprops=dict(arrowstyle='->', lw=0.7, color='k', shrinkA=2, shrinkB=2))
ax_d.text(0.97, 0.04, 'Arrow: more rain needed\nfor the same drying time', transform=ax_d.transAxes, ha='right', va='bottom', fontsize=6.5)
# ---- e: rain needed vs rain gain expected from warming
vals = [('To restore\n1979–2000\ndrying time', FIT.factor, FIT.f_lo, FIT.f_hi, OI['verm']), ('From\nwarming\n(7% per K)', 1.07 ** FIT.dT, None, None, '#7F7F7F')]
for i, (lb, v, lo, hi, c) in enumerate(vals):
    ax_e.bar(i, 100 * (v - 1), color=c, width=0.6, alpha=0.9)
    if lo: ax_e.plot([i, i], [100 * (lo - 1), 100 * (hi - 1)], color='k', lw=0.7)
    ax_e.text(i, 100 * ((hi or v) - 1) + 3, f'+{100 * (v - 1):.0f}%', ha='center', va='bottom', fontsize=6.5, color=c)
ax_e.set_xticks([0, 1]); ax_e.set_xticklabels([v[0] for v in vals], fontsize=6.5); ax_e.set_xlim(-0.6, 1.6); ax_e.set_ylim(0, 115)
ax_e.set_ylabel('Increase in rain amount (%)', labelpad=1); ax_e.tick_params(axis='x', length=0)
for s, xx, yy in [('a', 0.5, HMM - 0.5), ('b', 120, HMM - 0.5), ('c', 0.5, my - 12.5), ('d', 60, my - 12.5), ('e', 133, my - 12.5)]:
    fig.text(xx / WMM, yy / HMM, s, fontsize=9, fontweight='bold', va='top', ha='left')
bad = overlap_check(fig); print('overlap check:', bad or 'none')
fig.savefig(HERE + '/Fig3.png', dpi=600); fig.savefig(HERE + '/Fig3.pdf')
from PIL import Image; im = Image.open(HERE + '/Fig3.png'); im.resize((im.size[0] // 3, im.size[1] // 3)).save(HERE + '/Fig3_prev.png')
