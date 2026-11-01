# Shared data and helpers for the final figure set (rain-bought fire-safe days).
import numpy as np, pandas as pd, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nature.mplstyle'))
# CEE: one sans typeface (Arial; Liberation Sans is metric-compatible), optimum 8 pt; ticks/annotations 7 pt, panel letters 9 pt bold
import matplotlib.font_manager as _fm
_SANS = 'Arial' if any(f.name == 'Arial' for f in _fm.fontManager.ttflist) else 'Liberation Sans'
plt.rcParams.update({'font.sans-serif': [_SANS, 'DejaVu Sans'], 'font.size': 7, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 7, 'axes.labelpad': 2, 'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5, 'legend.handletextpad': 0.4,
                     'legend.borderaxespad': 0.3, 'legend.labelspacing': 0.3, 'pdf.fonttype': 42})
FS = dict(tick=7, label=8, ann=7, letter=9)
def letter(ax, s, x, y, fig=None, **k):
    (fig or ax.figure).text(x, y, s, fontsize=9, fontweight='bold', va='top', ha='left', **k) if fig is not None else ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', ha='left', **k)
C = './'
MM = 1 / 25.4
# Okabe-Ito
OI = dict(orange='#E69F00', sky='#56B4E9', green='#009E73', yellow='#F0E442', blue='#0072B2', verm='#D55E00', purple='#CC79A7', grey='#999999', black='#000000')
STATION_REGS = ['Mediterranean', 'US West', 'Central Asia', 'Mexico', 'Australia', 'Siberia', 'Canada West']
GRID_REGS = ['S America', 'Sahel/C Africa', 'S Africa', 'SE Asia', 'India']
NAME = {'S America': 'South America', 'Sahel/C Africa': 'Sahel & Central Africa', 'S Africa': 'Southern Africa', 'SE Asia': 'Southeast Asia',
        'Mediterranean': 'Mediterranean', 'US West': 'Western US', 'Central Asia': 'Central Asia', 'Mexico': 'Mexico', 'Australia': 'Australia',
        'Siberia': 'Siberia', 'Canada West': 'Western Canada', 'India': 'India'}

import matplotlib.patheffects as _pe
HALO = [_pe.withStroke(linewidth=2.2, foreground='white')]

def binned(base, bounds):
    # discrete classes for change maps: labelled bin edges, arrow ends for values beyond the last bin
    cols = plt.get_cmap(base)(np.linspace(0.01, 0.99, len(bounds) + 1))
    cm = matplotlib.colors.ListedColormap(cols[1:-1]); cm.set_under(cols[0]); cm.set_over(cols[-1])
    return cm, matplotlib.colors.BoundaryNorm(bounds, len(bounds) - 1)

def rd(f):
    R = pd.read_csv(C + f); R['rel'] = 100 * (R.rec_l / R.rec_e - 1); R['dRH'] = R.H_l - R.H_e; R['dT'] = R.T_l - R.T_e
    if 'spn_e' in R: R['prot_e'] = R.spn_e * R.spl_e; R['prot_l'] = R.spn_l * R.spl_l
    return R

def stations(): return rd('ffmc5_station_era5_thr10.csv')
def grid(): return rd('ffmc7_station_era5_chirps_obs_grid_thr10.csv')

def block_ci(df, col='rel', stat=np.median, nb=2000, size=5., seed=0):
    """95% CI of a regional statistic by resampling 5-degree spatial blocks (accounts for spatial correlation)."""
    rng = np.random.default_rng(seed); b = (np.floor(df.lat / size).astype(int).astype(str) + '_' + np.floor(df.lon / size).astype(int).astype(str)).values
    groups = [df[col].values[b == u] for u in np.unique(b)]; out = []
    for _ in range(nb):
        pick = rng.integers(0, len(groups), len(groups)); out.append(stat(np.concatenate([groups[i] for i in pick])))
    return np.percentile(out, [2.5, 97.5])

def region_table():
    S, G = stations(), grid(); rows = []
    for R, regs, kind in [(S, STATION_REGS, 'station'), (G, GRID_REGS, 'grid')]:
        for r in regs:
            g = R[R.reg == r]
            if len(g) < 10: continue
            if NAME[r] in {'Mexico', 'Western Canada', 'Central Asia'}:   # few stations: resample stations, not 5-degree blocks
                rng = np.random.default_rng(1); v = g.rel.values
                bs = [np.median(rng.choice(v, len(v))) for _ in range(4000)]
                lo, hi = np.percentile(bs, [2.5, 97.5])
            else:
                lo, hi = block_ci(g)
            rows.append(dict(reg=r, name=NAME[r], kind=kind, n=len(g), rel=g.rel.median(), lo=lo, hi=hi, shorter=(g.rel < 0).mean(),
                             dRH=g.dRH.median(), dT=g.dT.median(), prot_e=g.prot_e.median(), prot_l=g.prot_l.median()))
    return pd.DataFrame(rows)

def trend_test(s):
    from scipy import stats
    s = s.dropna(); x = s.index.values.astype(float); y = s.values
    bt = np.polyfit(x, y, 1); res = y - np.polyval(bt, x)
    r1 = np.corrcoef(res[:-1], res[1:])[0, 1]; n = len(y); neff = max(n * (1 - r1) / (1 + r1), 3)
    se = np.sqrt(np.sum(res ** 2) / (neff - 2) / np.sum((x - x.mean()) ** 2))
    return bt[0] * 10, 2 * stats.t.sf(abs(bt[0] / se), neff - 2)

def label(ax, s, x=-0.02, y=1.0, fig=None):
    (fig or ax.figure).text(x, y, s, fontsize=8, fontweight='bold', va='bottom', ha='left', transform=ax.transAxes if fig is None else None)

def overlap_check(fig):
    """Report overlapping text boxes (axes labels, ticks, annotations) in a rendered figure."""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); boxes = []
    undrawn = set()
    for a in fig.axes:
        for axis in (a.xaxis, a.yaxis):
            drawn = set(id(k) for k in axis._update_ticks())
            for k in axis.get_major_ticks() + axis.get_minor_ticks():
                if id(k) not in drawn: undrawn.update([id(k.label1), id(k.label2)])
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text().strip() or id(t) in undrawn: continue
        try: bb = t.get_window_extent(r)
        except Exception: continue
        if bb.width < 1: continue
        boxes.append((t.get_text()[:30], bb, t.get_rotation() not in (0, 90, 270) and not t.get_text().startswith('\n')))
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if boxes[i][2] and boxes[j][2]: continue   # two rotated tick labels: their upright boxes overlap, the glyphs do not
            if a.overlaps(b) and min(a.x1, b.x1) - max(a.x0, b.x0) > 1.5 and min(a.y1, b.y1) - max(a.y0, b.y0) > 1.5: bad.append((boxes[i][0], boxes[j][0]))
    W, H = fig.bbox.width, fig.bbox.height
    out = [('OUTSIDE FIGURE', s) for s, bb, _ in boxes if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > W + 0.5 or bb.y1 > H + 0.5]
    return bad + out + [('OCCLUDED', h) for h in occlusion_check(fig)]

def occlusion_check(fig, pad=0.5):
    """Report text (annotations, legends, labels) crossed by data lines, error bars or scatter markers in any axes."""
    from matplotlib.collections import PathCollection, LineCollection
    from matplotlib.legend import Legend
    fig.canvas.draw(); r = fig.canvas.get_renderer(); dpi = fig.dpi
    boxes = []
    for leg in fig.findobj(Legend):
        if leg.get_visible(): boxes.append(('legend:' + ','.join(t.get_text()[:12] for t in leg.get_texts())[:40], leg.get_window_extent(r), leg))
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text().strip(): continue
        if any(isinstance(p, Legend) for p in [t.axes.get_legend() if t.axes else None] + fig.legends if p is not None and t in p.get_texts()): continue
        if getattr(t, '_is_tick', False): continue
        try: bb = t.get_window_extent(r)
        except Exception: continue
        if bb.width < 1: continue
        boxes.append((t.get_text()[:30], bb, t))
    def pts_line(ln):
        xy = ln.get_transform().transform(ln.get_xydata()); xy = xy[np.isfinite(xy).all(1)]
        if len(xy) < 2: return xy
        out = [xy[:1]]
        for a, b in zip(xy[:-1], xy[1:]):
            n = int(max(2, np.hypot(*(b - a)) / 1.0)); out.append(a + (b - a) * np.linspace(0, 1, n)[:, None])
        return np.vstack(out)
    hits = []
    for ax in fig.axes:
        if not ax.axison: art_lines = []; colls = [c for c in ax.collections if isinstance(c, PathCollection)]
        else: art_lines = [l for l in ax.lines if l.get_visible() and l.get_linestyle() not in ('None', '', ' ')]; colls = [c for c in ax.collections if isinstance(c, (PathCollection, LineCollection))]
        ctx = ax.get_window_extent(r)
        geoms = []
        for l in art_lines: geoms.append(('line', pts_line(l), 0.5 * l.get_linewidth() * dpi / 72))
        for l in ax.lines:
            if l.get_marker() not in (None, 'None', '', ' ') and l.get_visible():
                xy = l.get_transform().transform(l.get_xydata()); geoms.append(('marker', xy[np.isfinite(xy).all(1)], 0.5 * l.get_markersize() * dpi / 72))
        for c in colls:
            if isinstance(c, PathCollection):
                off = c.get_offsets()
                if len(off) == 0: continue
                xy = c.get_offset_transform().transform(off); s = c.get_sizes(); rad = (np.sqrt(s.max() if len(s) else 20) / 2) * dpi / 72
                geoms.append(('scatter', xy[np.isfinite(xy).all(1)], rad))
            else:
                for seg in c.get_segments():
                    if len(seg) < 2: continue
                    xy = c.get_transform().transform(seg); out = []
                    for a, b in zip(xy[:-1], xy[1:]):
                        n = int(max(2, np.hypot(*(b - a)) / 1.0)); out.append(a + (b - a) * np.linspace(0, 1, n)[:, None])
                    lw = (c.get_linewidths()[0] if len(c.get_linewidths()) else 1) * dpi / 72 / 2; geoms.append(('errbar', np.vstack(out), lw))
        for name, bb, art in boxes:
            if art in (ax.xaxis.label, ax.yaxis.label) or art in ax.get_xticklabels() + ax.get_yticklabels(): continue
            if not bb.overlaps(ctx): continue
            halo = any('Stroke' in type(pe).__name__ for pe in (art.get_path_effects() if hasattr(art, 'get_path_effects') else []))
            for kind, xy, rad in geoms:
                if len(xy) == 0 or (halo and kind in ('line', 'errbar')): continue
                inside = (xy[:, 0] > bb.x0 - rad + pad) & (xy[:, 0] < bb.x1 + rad - pad) & (xy[:, 1] > bb.y0 - rad + pad) & (xy[:, 1] < bb.y1 + rad - pad)
                if inside.any(): hits.append((name, kind)); break
    return sorted(set(hits))
