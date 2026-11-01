# Fig. 1 Fine fuels dry faster after rain across most fire-prone land (ERA5 1.5-degree grid, 1979-2000 vs 2001-2022).
import json, numpy as np, pandas as pd, sys, os, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import plt, MM, OI, HALO, binned, overlap_check, trend_test
import cartopy.crs as ccrs, cartopy.feature as cfeature
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/data/'
R = pd.read_csv(D + 'fig1_cells.csv'); G = pd.read_csv(D + 'fig1_regions.csv'); S = pd.read_csv(D + 'fig1_series.csv', index_col=0)
RB = pd.read_csv(D + 'fig1_robust.csv'); JS = json.load(open(D + 'figS5_summary.json')); JT = json.load(open(D + 'figS5_txrn.json'))
m = lambda v, f='.1f': format(v, f).replace('-', '−')
BOUNDS = [-30, -20, -15, -10, -5, -2, 2, 5, 10, 15, 20, 30]
CM, NORM = binned('BrBG', BOUNDS)
col = lambda v: CM(NORM(v)) if BOUNDS[0] <= v <= BOUNDS[-1] else (CM.get_under() if v < BOUNDS[0] else CM.get_over())
PC = ccrs.PlateCarree(); ROB = ccrs.Robinson(central_longitude=10)
LAT0, LAT1 = -57, 80
yproj = lambda lat: ROB.transform_points(PC, np.zeros_like(np.asarray(lat, float)), np.asarray(lat, float))[:, 1]
Y0, Y1 = yproj([LAT0, LAT1]); X0, X1 = ROB.x_limits
WMM, HMM = 183, 116
fig = plt.figure(figsize=(WMM * MM, HMM * MM))
ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)   # mm from lower-left
mw = 117; mh = mw * (Y1 - Y0) / (X1 - X0); my = HMM - 1.5 - mh
ax_map = ax_(1, my, mw, mh, projection=ROB)
ax_z = ax_(mw + 3, my, 10, mh)
ax_r = ax_(146, my - 3, 36, mh + 3.5)
ax_b = ax_(15, 9.5, 97, 37)
ax_c = ax_(129, 9.5, 53, 37)
# ---- a: map
ax_map.set_xlim(X0, X1); ax_map.set_ylim(Y0, Y1); ax_map.set_xticks([]); ax_map.set_yticks([])
ax_map.add_feature(cfeature.LAND.with_scale('110m'), facecolor='#D9D9D9', edgecolor='none', zorder=0)
grid = np.full((121, 240), np.nan)
ia = np.rint((R.lat.values + 90) / 1.5).astype(int); io = np.rint((R.lon.values % 360) / 1.5).astype(int) % 240
grid[ia, io] = R.rel.values
lon_e = np.arange(241) * 1.5 - 0.75; lat_e = np.clip(np.arange(122) * 1.5 - 90.75, -90, 90)
ax_map.pcolormesh(lon_e, lat_e, grid, cmap=CM, norm=NORM, transform=PC, zorder=1, rasterized=True)
ax_map.coastlines('110m', linewidth=0.3, color='#4D4D4D', zorder=2)
ax_map.spines['geo'].set_linewidth(0.5)
ax_map.gridlines(crs=PC, linewidth=0.3, color='#B3B3B3', linestyle='--', xlocs=range(-180, 181, 60), ylocs=[-30, 0, 30, 60], zorder=1.5)
# inset histogram of the cell changes (map colours), in the South Pacific
ih = ax_(1 + mw * 0.05 + 3.5, my + mh * 0.12 - 1, mw * 0.155, mh * 0.25)
v = np.clip(R.rel.values, -45, 45); edges = np.arange(-45, 46, 5)
hc, _ = np.histogram(v, edges); mids = 0.5 * (edges[:-1] + edges[1:])
ih.bar(mids, hc, width=4.4, color=[col(q) for q in mids], edgecolor='none')
ih.axvline(0, color='k', lw=0.5); ih.set_xlim(-45, 45); ih.set_yticks([]); ih.spines['left'].set_visible(False)
ih.set_xticks([-40, 0, 40]); ih.set_xticklabels(['−40', '0', '40']); ih.tick_params(labelsize=6, length=2, pad=1); ih.patch.set_alpha(0)
ih.text(0.0, 1.30, f'Decline in\n{(R.rel < 0).mean() * 100:.0f}% of cells', transform=ih.transAxes, fontsize=6.5, color=col(-17), ha='left', va='bottom', path_effects=HALO)
ih.text(0.0, 1.04, f'Median {m(R.rel.median())}%', transform=ih.transAxes, fontsize=6.5, ha='left', va='bottom', path_effects=HALO)
# colour bar under the map
cax = ax_(22, my - 5.2, 76, 2.2)
cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=NORM, cmap=CM), cax=cax, orientation='horizontal', extend='both', ticks=BOUNDS)
cb.ax.set_xticklabels([m(b, 'd') for b in BOUNDS]); cb.ax.tick_params(labelsize=6.5, length=2, pad=1); cb.outline.set_linewidth(0.4)
cb.set_label('Change in drying time after rain, 2001–2022 vs 1979–2000 (%)', fontsize=7, labelpad=1.5)
# ---- zonal median strip in the map's Robinson y coordinates
z = []
for la in np.sort(R.lat.unique()):
    s = R[np.abs(R.lat - la) <= 2.3].rel
    if len(s) >= 8: z.append((la, s.median(), s.quantile(0.25), s.quantile(0.75)))
z = np.array(z); yz = yproj(z[:, 0])
ax_z.set_ylim(Y0, Y1)
ax_z.fill_betweenx(yz, z[:, 2], z[:, 3], color='#BDBDBD', lw=0, alpha=0.7)
ax_z.plot(z[:, 1], yz, color='k', lw=0.8); ax_z.axvline(0, color='k', lw=0.4, ls='--')
ax_z.set_yticks(yproj([-30., 0., 30., 60.])); ax_z.set_yticklabels(['30° S', '0°', '30° N', '60° N'])
ax_z.yaxis.tick_right(); ax_z.spines['left'].set_visible(False); ax_z.spines['right'].set_visible(True)
ax_z.set_xlim(-50, 15); ax_z.set_xticks([-40, 0]); ax_z.set_xticklabels(['−40', '0'])
ax_z.tick_params(labelsize=6.5, pad=1)
ax_z.set_xlabel('Zonal (%)', fontsize=7, labelpad=1)
# ---- b: robustness (labels above each interval)
RB = RB.iloc[::-1].reset_index(drop=True)
SHORT = {'Reference (ERA5 weather and rain)': 'Reference (ERA5)', 'Satellite–gauge rain (CHIRPS)': 'CHIRPS rain', 'Gauge rain (CPC)': 'CPC gauge rain',
         'Observed noon wind': 'Daily ERA5 noon wind', 'Daily maximum T, minimum RH': 'Daily max. T, min. RH', 'Dry threshold FFMC 90': 'Threshold FFMC 90',
         'Nelson model, 1-h stick': 'Nelson model, 1-h fuel', 'Nelson model, 10-h stick': 'Nelson model, 10-h fuel',
         'Weather stations, station rain': '510 stations, station rain', 'ISD stations, measured humidity': '555 stations, measured RH',
         'ISD stations, measured humidity, CPC rain': 'Measured RH and CPC rain',
         'ISD stations, station only, daily max/min': '300 stations, station only'}
grp = {'Reference': 'k', 'Satellite': OI['blue'], 'Gauge': OI['blue'], 'Nelson': OI['green'], 'Weather': OI['purple'], 'ISD': OI['purple']}
ax_r.axvline(0, color='k', lw=0.5); ax_r.axvline(RB.med.iloc[-1], color='#9E9E9E', lw=0.5, ls='--')
for i, r in RB.iterrows():
    i = i + 0.4 * min(i, 3)                       # extra space for the three station rows with an ERA5 marker
    c = next((v for k, v in grp.items() if r.version.startswith(k)), '#6E6E6E')
    ax_r.plot([r.lo, r.hi], [i, i], color=c, lw=0.9, solid_capstyle='butt'); ax_r.scatter(r.med, i, s=12, color=c, zorder=3, edgecolor='none')
    ax_r.text(-27.3, i + 0.22, SHORT[r.version], fontsize=6.5, color=c, ha='left', va='bottom', path_effects=HALO)
    if r.version.startswith('ISD'):                # ERA5 temperature and humidity at the same stations and hours (open marker)
        q = JT['B_all']['era'] if 'station only' in r.version else JS['cpc' if 'CPC' in r.version else 'main']['all']['era']
        ax_r.plot([q[1], q[2]], [i - 0.3, i - 0.3], color='#8C8C8C', lw=0.7, solid_capstyle='butt')
        ax_r.scatter(q[0], i - 0.3, s=11, facecolor='white', edgecolor=c, linewidth=0.7, zorder=3)
ax_r.set_ylim(-0.55, len(RB) - 1 + 1.2 + 0.85); ax_r.set_yticks([]); ax_r.spines['left'].set_visible(False)
ax_r.set_xlim(-28, 2); ax_r.set_xticks([-20, -10, 0]); ax_r.set_xticklabels(['−20', '−10', '0'])
ax_r.set_xlabel('Median change (%)', labelpad=1)
# ---- c: AR6 regions
G = G.reset_index(drop=True); x = np.arange(len(G))
ax_b.axhline(0, color='k', lw=0.5); ax_b.axhline(R.rel.median(), color='#9E9E9E', lw=0.6, ls='--')
for i, r in G.iterrows():
    ax_b.plot([i, i], [r.lo, r.hi], color='#595959', lw=0.7, solid_capstyle='butt', zorder=1)
    ax_b.scatter(i, r.med, s=16, color=col(r.med), edgecolor='k', linewidth=0.4, zorder=2)
ax_b.set_xticks(x); ax_b.set_xticklabels(G.ar6, rotation=90, fontsize=6.5); ax_b.set_xlim(-0.7, len(G) - 0.3)
ax_b.set_ylim(-45, 20); ax_b.set_yticks([-40, -20, 0, 20]); ax_b.set_yticklabels(['−40', '−20', '0', '20'])
ax_b.set_ylabel('Change in drying\ntime after rain (%)', labelpad=2)
ax_b.text(len(G) - 0.6, -41, f'All cells {m(R.rel.median())}% (dashed line)', ha='right', va='bottom', fontsize=6.5, color='#595959')
ax_b.tick_params(axis='x', pad=1)
# ---- d: global series
for k, c, lab, yy in [('era5', 'k', 'ERA5 rain', 0.05), ('chirps', OI['blue'], 'CHIRPS rain', 0.16)]:
    s = S[k].dropna(); ax_c.plot(s.index, s.values, color=c, lw=0.9)
    b = np.polyfit(s.index, s.values, 1); ax_c.plot(s.index, np.polyval(b, s.index), color=c, lw=0.6, ls='--')
    tr, p = trend_test(s)
    ax_c.text(0.03, yy, f'{lab}: {m(tr)}% per decade', transform=ax_c.transAxes, color=c, fontsize=6.5, ha='left', va='bottom')
ax_c.axhline(0, color='k', lw=0.4)
ax_c.set_xlim(1978, 2023); ax_c.set_xticks([1980, 1990, 2000, 2010, 2020]); ax_c.set_ylim(-16, 16)
ax_c.set_yticks([-10, 0, 10]); ax_c.set_yticklabels(['−10', '0', '10'])
ax_c.set_ylabel('Anomaly of drying\ntime after rain (%)', labelpad=2)
ax_c.text(0.97, 0.97, '$P$ < 0.001 (both)', transform=ax_c.transAxes, ha='right', va='top', fontsize=6.5)
# ---- panel letters (shared height per row)
for s, xx, yy in [('a', 0.5, HMM - 0.5), ('b', 141, HMM - 0.5), ('c', 0.5, 50.5), ('d', 116, 50.5)]:
    fig.text(xx / WMM, yy / HMM, s, fontsize=9, fontweight='bold', va='top', ha='left')
bad = overlap_check(fig); print('overlap check:', bad or 'none')
fig.savefig(HERE + '/Fig1.png', dpi=600); fig.savefig(HERE + '/Fig1.pdf')
from PIL import Image; im = Image.open(HERE + '/Fig1.png'); im.resize((im.size[0] // 3, im.size[1] // 3)).save(HERE + '/Fig1_prev.png')
