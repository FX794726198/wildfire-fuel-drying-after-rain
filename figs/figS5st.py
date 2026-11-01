# Fig. S5 Station humidity check: drying time after rain from measured noon temperature and humidity at 555 ISD stations,
# compared with ERA5 at the same cell and hour (same rain and fire-season months).
import numpy as np, pandas as pd, sys, os, json, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import plt, MM, OI, HALO, binned, overlap_check
from figmaps import ROB, PC, ASPECT, X0, X1, Y0, Y1, hcolorbar
import cartopy.feature as cfeature
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/data/'
m = lambda v, f='.1f': format(v, f).replace('-', '−')
X = pd.read_csv(D + 'figS5_stations_main.csv'); J = json.load(open(D + 'figS5_summary.json')); JT = json.load(open(D + 'figS5_txrn.json'))
BOUNDS = [-30, -20, -15, -10, -5, -2, 2, 5, 10, 15, 20, 30]; CM, NORM = binned('BrBG', BOUNDS)
WMM = 183; mw = 112; mh = mw * ASPECT; HMM = 134
fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
my = HMM - 1 - mh
a = ax_(1, my, mw, mh, projection=ROB)
b = ax_(132, my + 1, 49, mh - 1)
c = ax_(13, 8, 46, 50)
d = ax_(110, 8, 71, 58)
# ---- a: station map
a.set_xlim(X0, X1); a.set_ylim(Y0, Y1); a.set_xticks([]); a.set_yticks([])
a.add_feature(cfeature.LAND.with_scale('110m'), facecolor='#D9D9D9', edgecolor='none', zorder=0)
a.coastlines('110m', linewidth=0.3, color='#4D4D4D', zorder=2); a.spines['geo'].set_linewidth(0.5)
a.gridlines(crs=PC, linewidth=0.3, color='#B3B3B3', linestyle='--', xlocs=range(-180, 181, 60), ylocs=[-30, 0, 30, 60], zorder=1.5)
o = np.argsort(np.abs(X.st.values))
a.scatter(X.lon.values[o], X.lat.values[o], c=X.st.values[o], cmap=CM, norm=NORM, s=5, edgecolor='k', linewidth=0.15, transform=PC, zorder=3)
A = J['main']['all']
a.text(0.04, 0.30, f"{A['n']} stations, {A['countries']} countries\nMedian {m(A['st'][0])}%\nShorter at {A['shorter_st'] * 100:.0f}%",
       transform=a.transAxes, fontsize=6.5, ha='left', va='top', path_effects=HALO)
hcolorbar(fig, ax_(19, my - 5.2, 76, 2.2), CM, NORM, BOUNDS, 'Change in drying time after rain, station humidity (%)')
# ---- b: humidity change, station vs ERA5
lim = (-10, 8)
b.plot(lim, lim, color='#9E9E9E', lw=0.6, ls='--', zorder=1)
def edge(ax, x, y, lim):                           # stations beyond the axes drawn at the edge as open circles
    out = (x < lim[0]) | (x > lim[1]) | (y < lim[0]) | (y > lim[1]); k = 0.985 * np.array(lim)
    ax.scatter(x[~out], y[~out], s=4, color='#6E6E6E', alpha=0.6, linewidths=0, zorder=2)
    ax.scatter(np.clip(x[out], *k), np.clip(y[out], *k), s=7, facecolor='none', edgecolor='#6E6E6E', linewidth=0.5, zorder=2, clip_on=False)
    return int(out.sum())
nb = edge(b, X.dRH_era.values, X.dRH_st.values, lim)
b.scatter([A['dRH_era'][0]], [A['dRH_st'][0]], s=40, marker='D', color=OI['verm'], edgecolor='k', linewidth=0.4, zorder=3)
b.set_xlim(*lim); b.set_ylim(*lim); b.set_xticks([-8, -4, 0, 4]); b.set_yticks([-8, -4, 0, 4])
b.set_xticklabels(['−8', '−4', '0', '4']); b.set_yticklabels(['−8', '−4', '0', '4'])
b.set_xlabel('RH change, ERA5 same hour (%)', labelpad=1); b.set_ylabel('RH change, station (%)', labelpad=1)
b.text(0.97, 0.03, f"Median: station {m(A['dRH_st'][0])},\nERA5 {m(A['dRH_era'][0])}\n$r$ = {m(A['r_dRH'], '.2f')}", transform=b.transAxes, fontsize=6.5, ha='right', va='bottom')
# ---- c: drying-time change, station vs ERA5
lim = (-60, 60)
c.plot(lim, lim, color='#9E9E9E', lw=0.6, ls='--', zorder=1); c.axhline(0, color='k', lw=0.3); c.axvline(0, color='k', lw=0.3)
nc = edge(c, X.era.values, X.st.values, lim); print('beyond axes: b', nb, 'c', nc)
c.scatter([A['era'][0]], [A['st'][0]], s=40, marker='D', color=OI['verm'], edgecolor='k', linewidth=0.4, zorder=3)
c.set_xlim(*lim); c.set_ylim(*lim); c.set_xticks([-50, 0, 50]); c.set_yticks([-50, 0, 50])
c.set_xticklabels(['−50', '0', '50']); c.set_yticklabels(['−50', '0', '50'])
c.set_xlabel('Drying-time change, ERA5 (%)', labelpad=1); c.set_ylabel('Drying-time change,\nstation (%)', labelpad=1)
c.text(0.97, 0.03, f"$r$ = {m(A['r_d'], '.2f')}", transform=c.transAxes, fontsize=6.5, ha='right', va='bottom')
# ---- d: medians with 95% intervals, station (filled) vs ERA5 at the same stations (open)
rows = [('hdr', None, 'Noon temperature and humidity'),
        ('main', 'all', 'Gaps filled with ERA5, ERA5 rain'), ('cpc', 'all', 'Gaps filled with ERA5, CPC rain'), ('main', 'cov90', '≥90% days measured'),
        ('main', 'no_break', 'No break in RH difference'), ('main', 'non_US', 'Outside the United States'),
        ('main', 'N extratropics', 'Northern extratropics'), ('main', 'S extratropics', 'Southern extratropics'), ('main', 'tropics', 'Tropics'),
        ('hdr', None, 'Daily max. temperature, min. humidity'),
        ('txrn', 'A_all', 'Gaps filled with ERA5'), ('txrn', 'B_all', 'Station only'), ('txrn', 'B_non_US', 'Station only, outside the US'),
        ('txrn', 'B_N extratropics', 'Station only, N extratropics'), ('txrn', 'B_S extratropics', 'Station only, S extratropics'), ('txrn', 'B_tropics', 'Station only, tropics')]
d.axvline(0, color='k', lw=0.5); XL = (-20, 15)
for i, (t, g, lab) in enumerate(rows[::-1]):
    if t == 'hdr':
        d.text(-0.02, i, lab, fontsize=6.5, fontweight='bold', ha='right', va='center', transform=d.get_yaxis_transform()); continue
    q = JT[g] if t == 'txrn' else J[t][g]
    for key, dy, fc in [('st', 0.13, OI['purple']), ('era', -0.13, 'white')]:
        v = q[key]; ec = OI['purple'] if key == 'st' else '#595959'
        d.plot([v[1], v[2]], [i + dy, i + dy], color=ec, lw=0.8, solid_capstyle='butt')
        if v[0] > XL[1]:
            d.scatter(XL[1] - 0.3, i + dy, s=13, marker='>', color=fc, edgecolor=ec, linewidth=0.6, zorder=3, clip_on=False)
            d.text(XL[1] - 1.2, i + dy + 0.1, m(v[0]), fontsize=6, ha='right', va='bottom', color=ec)
        else: d.scatter(v[0], i + dy, s=13, marker='o', color=fc, edgecolor=ec, linewidth=0.6, zorder=3)
    d.text(-0.02, i, f"{lab} ({q['n']})", fontsize=6.5, ha='right', va='center', transform=d.get_yaxis_transform())
d.set_ylim(-0.6, len(rows) - 0.4); d.set_yticks([]); d.spines['left'].set_visible(False)
d.set_xlim(*XL); d.set_xticks([-20, -10, 0, 10]); d.set_xticklabels(['−20', '−10', '0', '10'])
d.set_xlabel('Median change in drying time (%)', labelpad=1)
d.scatter([], [], s=13, color=OI['purple'], label='Station weather'); d.scatter([], [], s=13, color='white', edgecolor='#595959', linewidth=0.6, label='ERA5, same cells and hours')
d.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=6.5, handletextpad=0.2, borderaxespad=0.1, columnspacing=1.2)
for s, xx, yy in [('a', 0.5, HMM - 0.5), ('b', 120, HMM - 0.5), ('c', 0.5, 73), ('d', 63, 73)]:
    fig.text(xx / WMM, yy / HMM, s, fontsize=9, fontweight='bold', va='top', ha='left')
bad = overlap_check(fig); print('S5 overlap check:', bad or 'none')
fig.savefig(HERE + '/FigS5.png', dpi=600); fig.savefig(HERE + '/FigS5.pdf')
from PIL import Image; im = Image.open(HERE + '/FigS5.png'); im.resize((im.size[0] // 3, im.size[1] // 3)).save(HERE + '/FigS5_prev.png')
