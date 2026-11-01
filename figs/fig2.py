# Fig. 2 Satellite fire detections follow the time fuels take to dry after rain (MODIS 2001-2022, ERA5 FFMC).
import numpy as np, pandas as pd, sys, os, json, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import plt, MM, OI, HALO, binned, overlap_check
import cartopy.crs as ccrs, cartopy.feature as cfeature
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/data/'; tag = 'modis'
A = pd.read_csv(D + f'fig2_cells_{tag}.csv'); G = pd.read_csv(D + f'fig2_regions_{tag}.csv'); DO = pd.read_csv(D + f'fig2_dose_{tag}.csv')
S = json.load(open(f'./res/fire_summary_{tag}.json')); CI = json.load(open(D + f'fig2_ci_{tag}.json'))
m = lambda v, f='.2f': format(v, f).replace('-', '−')
BOUNDS = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5]
_seq = plt.get_cmap('YlGnBu')                                                     # sequential below 1 (fewer fires on wet days), warm hue above 1
CM = matplotlib.colors.ListedColormap(list(_seq(np.linspace(0.92, 0.14, 7))) + ['#F4A582']); CM.set_under(_seq(1.0)); CM.set_over('#D6604D')
NORM = matplotlib.colors.BoundaryNorm(BOUNDS, len(BOUNDS) - 1)
col = lambda v: CM(NORM(v)) if BOUNDS[0] <= v <= BOUNDS[-1] else (CM.get_under() if v < BOUNDS[0] else CM.get_over())
PC = ccrs.PlateCarree(); ROB = ccrs.Robinson(central_longitude=10)
yproj = lambda lat: ROB.transform_points(PC, np.zeros_like(np.asarray(lat, float)), np.asarray(lat, float))[:, 1]
Y0, Y1 = yproj([-57, 80]); X0, X1 = ROB.x_limits
WMM, HMM = 183, 146; R3 = 42                                               # R3: height of the new bottom row (cloud-corrected test)
fig = plt.figure(figsize=(WMM * MM, HMM * MM))
ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
mw = 117; mh = mw * (Y1 - Y0) / (X1 - X0); my = HMM - 1.5 - mh
ax_map = ax_(1, my, mw, mh, projection=ROB)
ax_d = ax_(133, my + 1, 49, mh - 12); axn = ax_(133, my + mh - 9, 49, 8.5)
ax_c = ax_(15, R3 + 9.5, 97, my - 22 - R3)
ax_e = ax_(151, R3 + 9.5, 31, my - 22 - R3)
# ---- a: map of wet/dry fire-rate ratios
ax_map.set_xlim(X0, X1); ax_map.set_ylim(Y0, Y1); ax_map.set_xticks([]); ax_map.set_yticks([])
ax_map.add_feature(cfeature.LAND.with_scale('110m'), facecolor='#D9D9D9', edgecolor='none', zorder=0)
grid = np.full((121, 240), np.nan)
ia = np.rint((A.lat.values + 90) / 1.5).astype(int); io = np.rint((A.lon.values % 360) / 1.5).astype(int) % 240
grid[ia, io] = (A.num / A.den).values
lon_e = np.arange(241) * 1.5 - 0.75; lat_e = np.clip(np.arange(122) * 1.5 - 90.75, -90, 90)
ax_map.pcolormesh(lon_e, lat_e, grid, cmap=CM, norm=NORM, transform=PC, zorder=1, rasterized=True)
ax_map.coastlines('110m', linewidth=0.3, color='#4D4D4D', zorder=2); ax_map.spines['geo'].set_linewidth(0.5)
ax_map.gridlines(crs=PC, linewidth=0.3, color='#B3B3B3', linestyle='--', xlocs=range(-180, 181, 60), ylocs=[-30, 0, 30, 60], zorder=1.5)
rr = (A.num / A.den)
ax_map.text(0.075, 0.36, f'All cells: {m(S["MH_RR_all"])}\n{(rr < 1).mean() * 100:.0f}% of cells below 1', transform=ax_map.transAxes, fontsize=6.5, ha='left', va='top', path_effects=HALO)
cax = ax_(22, my - 5.2, 76, 2.2)
cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=NORM, cmap=CM), cax=cax, orientation='horizontal', extend='both', ticks=BOUNDS)
cb.ax.set_xticklabels([f'{b:g}' for b in BOUNDS]); cb.ax.tick_params(labelsize=6.5, length=2, pad=1); cb.outline.set_linewidth(0.4)
cb.set_label('Fire detections on days with wet fuels relative to days with dry fuels (rate ratio)', fontsize=7, labelpad=1.5)
# ---- b: dose-response
ax_d.fill_between(DO.rc, DO.lo, DO.hi, color=OI['blue'], alpha=0.2, lw=0)
ax_d.plot(DO.rc, DO.ratio, color=OI['blue'], lw=1.0, marker='o', ms=2.5)
ax_d.axhline(1, color='k', lw=0.5, ls='--')
ax_d.set_xlim(0.5, 15.5); ax_d.set_xticks([1, 5, 10, 15]); ax_d.set_xticklabels(['1', '5', '10', '≥15'])
ax_d.set_ylim(0, 1.1); ax_d.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_d.set_xlabel('Drying time after the rain (days)', labelpad=1); ax_d.set_ylabel('Fires in 10 days after rain\nrelative to matched dry spells', labelpad=2)
ax_d.text(0.97, 0.80, f'Each day of drying time: fewer fires,\n{m(S["slope"])} ({m(S["slope_ci"][0])}–{m(S["slope_ci"][1])}) days of\naverage fire activity',
          transform=ax_d.transAxes, ha='right', va='top', fontsize=6.5, color=OI['blue'])
axn.bar(DO.rc, DO.n / 1000, color='#BDBDBD', width=0.75); axn.set_xlim(0.5, 15.5); axn.set_xticks([1, 5, 10, 15]); axn.set_xticklabels([])
axn.set_yticks([0, 10]); axn.tick_params(labelsize=6.5, length=2, pad=1)
axn.text(0.97, 0.92, 'Rains (thousands)', transform=axn.transAxes, fontsize=6.5, va='top', ha='right', color='#6E6E6E')
# ---- c: AR6 regions
G = G.reset_index(drop=True); x = np.arange(len(G))
ax_c.axhline(1, color='k', lw=0.5, ls='--')
for i, r in G.iterrows():
    ax_c.plot([i, i], [r.lo, r.hi], color='#595959', lw=0.7, solid_capstyle='butt')
    ax_c.scatter(i, r.rr, s=16, color=col(r.rr), edgecolor='k', linewidth=0.4, zorder=2)
ax_c.set_yscale('log'); ax_c.set_ylim(0.05, 1.5); ax_c.set_yticks([0.05, 0.1, 0.2, 0.5, 1]); ax_c.set_yticklabels(['0.05', '0.1', '0.2', '0.5', '1'])
ax_c.minorticks_off()
ax_c.set_xticks(x); ax_c.set_xticklabels(G.ar6, rotation=90, fontsize=6.5); ax_c.set_xlim(-0.7, len(G) - 0.3); ax_c.tick_params(axis='x', pad=1)
ax_c.set_ylabel('Wet/dry rate ratio', labelpad=2)
# ---- d: checks
CH = [('Wet vs dry days, same\ncell, year and month', S['MH_RR_all'], CI['MH_ci'], 'k'),
      ('Median over cells', S['MH_RR_cell_median'], CI['MH_cell_median_ci'], '#6E6E6E'),
      ('10 days after rain vs\nmatched dry spells', S['avoid_ratio'], CI['avoid_ratio_ci'], OI['blue'])]
for i, (lab, v, ci, c) in enumerate(CH[::-1]):
    ax_e.barh(i, v, color=c, height=0.55, alpha=0.85); ax_e.plot(ci, [i, i], color='k', lw=0.7, solid_capstyle='butt')
    if v > 0.5: ax_e.text(ci[0] - 0.03, i, m(v), va='center', ha='right', fontsize=6.5, color='white')
    else: ax_e.text(ci[1] + 0.02, i, m(v), va='center', ha='left', fontsize=6.5, color=c)
ax_e.set_yticks(range(len(CH))); ax_e.set_yticklabels([c[0] for c in CH[::-1]], fontsize=6.5)
ax_e.axvline(1, color='k', lw=0.5, ls='--'); ax_e.set_xlim(0, 1.1); ax_e.set_xticks([0, 0.5, 1])
ax_e.set_xlabel('Ratio of fire detections', labelpad=1); ax_e.tick_params(axis='y', length=0, pad=2)
# ---- e-g: fire activity on cloud-free land after rain, five regions together (e), by region (f) and by condition (g); Fig. S3 has the checks
FR = json.load(open(D + 'cloudfree_results.json')); get = lambda t: next(o for o in FR if o['tag'] == t)
CC = [OI['orange'], OI['sky'], OI['blue']]; T = np.arange(-10, 21)
def sm(v): v = np.asarray(v, float); v = np.where(np.isfinite(v), v, np.nanmean(v)); return np.convolve(v, np.ones(3) / 3, 'same')
ax = ax_(15, 8.5, 50, R3 - 10); o = get('five regions: main'); A_ = np.array(o['A'], float)
ax.axhline(1, color='#4D4D4D', lw=0.5, ls='--'); ax.axhline(0.5, color='#4D4D4D', lw=0.5, ls=':'); ax.axvline(0, ymax=0.6, color='#4D4D4D', lw=0.5)
jj = (T >= -7) & (T <= 15)
for cl in range(3):
    ax.plot(T[jj], sm(A_[cl])[jj], color=CC[cl], lw=1.0, label=['1 day', '2 days', '3–15 days'][cl])
    ax.plot(o['T50'][cl], 0.5, '|', ms=5, mew=1.2, color=CC[cl])
ax.set_xlim(-7, 15); ax.set_xticks([-5, 0, 5, 10, 15]); ax.set_xticklabels(['−5', '0', '5', '10', '15'])
ax.set_ylim(0, 1.75); ax.set_yticks([0, 0.5, 1, 1.5]); ax.set_xlabel('Days from rain', labelpad=1)
ax.set_ylabel('Fire activity relative to\npre-rain dry days (cloud-free)', labelpad=2)
ax.text(0.03, 0.97, f'Five regions, 2003–2022\nΔT50 {m(o["dT50"], ".1f")} days ({m(o["dT50_ci"][0], ".1f")} to {m(o["dT50_ci"][1], ".1f")})',
        transform=ax.transAxes, fontsize=6.5, ha='left', va='top', linespacing=1.3)
ax.legend(loc='upper right', frameon=False, title='Drying time', title_fontsize=6.5, fontsize=6.5, handlelength=1.2, borderaxespad=0.1, labelspacing=0.2)
# f: by region; left ΔT50, right fire activity on days 8-15 (all classes)
REG = [('five regions: main', 'All five'), ('region SAM', 'SAM'), ('region SES', 'SES'), ('region NAU', 'NAU'), ('region ESAF', 'ESAF'), ('region CAF', 'CAF')]
REGC = {'SAM': OI['green'], 'SES': OI['blue'], 'ESAF': OI['orange'], 'NAU': OI['purple'], 'CAF': OI['verm'], 'All five': 'k'}
f1 = ax_(84, 8.5, 29, R3 - 10); f2 = ax_(116, 8.5, 14, R3 - 10)
for i, (t, lb) in enumerate(REG):
    o = get(t); yy = len(REG) - 1 - i; c = REGC[lb]; cens = max(o['T50']) >= 21
    f1.plot(o['dT50_ci'], [yy, yy], color=c, lw=0.8, solid_capstyle='butt')
    f1.plot(o['dT50'], yy, 'o' if not cens else '>', ms=3.5 if i else 4.2, mfc=c if not cens else 'white', mec=c, mew=0.8)
    f2.plot(o['level_8_15'], yy, 'o', ms=3.5 if i else 4.2, color=c)
for a_ in (f1, f2):
    a_.set_ylim(-0.6, len(REG) - 0.4); a_.set_yticks(range(len(REG) - 1, -1, -1)); a_.tick_params(axis='y', length=0)
    a_.axhline(len(REG) - 1.5, color='#BFBFBF', lw=0.5)
f1.set_yticklabels([lb for _, lb in REG], fontsize=6.5)
for tl, (_, lb) in zip(f1.get_yticklabels(), REG): tl.set_color(REGC[lb])
f2.set_yticklabels([])
f1.axvline(0, color='k', lw=0.5, ls='--'); f1.set_xlim(-17, 21); f1.set_xticks([-10, 0, 10, 20]); f1.set_xticklabels(['−10', '0', '10', '20'])
f1.set_xlabel('ΔT50 (days)', labelpad=1)
f2.axvline(1, color='k', lw=0.5, ls='--'); f2.set_xlim(0.2, 1.45); f2.set_xticks([0.5, 1]); f2.set_xticklabels(['0.5', '1'])
f2.set_xlabel('Days 8–15', labelpad=1)
# g: by condition (five regions together)
from matplotlib.lines import Line2D
CND = [('season: first two-thirds', 'Early in season'), ('season: last third', 'Late in season'), (None, None)]
lr = [o['tag'] for o in FR if o['tag'].startswith('later rain')]; tc = [o['tag'] for o in FR if o['tag'].startswith('tree cover')]
CND += [(lr[0], 'Less later rain'), (lr[1], 'More later rain'), (None, None), (tc[0], 'Less tree cover'), (tc[1], 'More tree cover')]
g1 = ax_(158, 8.5, 24, R3 - 10)
for i, (t, lb) in enumerate(CND):
    if t is None: continue
    o = get(t); yy = len(CND) - 1 - i; cens = max(o['T50']) >= 21
    g1.plot(o['dT50_ci'], [yy, yy], color='k', lw=0.8, solid_capstyle='butt')
    g1.plot(o['dT50'], yy, 'o' if not cens else '>', ms=3.5, mfc='k' if not cens else 'white', mec='k', mew=0.8)
g1.set_ylim(-0.6, len(CND) - 0.4); g1.set_yticks([len(CND) - 1 - i for i, (t, _) in enumerate(CND) if t]); g1.tick_params(axis='y', length=0)
g1.set_yticklabels([lb for t, lb in CND if t], fontsize=6.5)
g1.axvline(0, color='k', lw=0.5, ls='--'); g1.set_xlim(-10, 21); g1.set_xticks([-10, 0, 10, 20]); g1.set_xticklabels(['−10', '0', '10', '20']); g1.set_xlabel('ΔT50 (days)', labelpad=1)
for s, xx, yy in [('a', 0.5, HMM - 0.5), ('b', 119, HMM - 0.5), ('c', 0.5, my - 9.5), ('d', 114, my - 9.5), ('e', 0.5, R3 + 0.5), ('f', 70, R3 + 0.5), ('g', 134, R3 + 0.5)]:
    fig.text(xx / WMM, yy / HMM, s, fontsize=9, fontweight='bold', va='top', ha='left')
bad = overlap_check(fig); print('overlap check:', bad or 'none')
fig.savefig(HERE + '/Fig2.png', dpi=600); fig.savefig(HERE + '/Fig2.pdf')
from PIL import Image; im = Image.open(HERE + '/Fig2.png'); im.resize((im.size[0] // 3, im.size[1] // 3)).save(HERE + '/Fig2_prev.png')
