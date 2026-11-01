# Supplementary figures S1-S4 and S6 (usage: python figS.py [S1 S2 ...]); S5 is figS5st.py, S7-S8 figS67.py.
import numpy as np, pandas as pd, sys, os, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import plt, MM, OI, HALO, binned, overlap_check
from figmaps import ROB, ASPECT, cellmap, hcolorbar
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/data/'
m = lambda v, f='.1f': format(v, f).replace('-', '−')
WANT = sys.argv[1:] or ['S1', 'S2', 'S3', 'S4', 'S6']
def save(fig, name):
    bad = overlap_check(fig); print(name, 'overlap check:', bad or 'none')
    fig.savefig(f'{HERE}/Fig{name}.png', dpi=600); fig.savefig(f'{HERE}/Fig{name}.pdf')
    from PIL import Image; im = Image.open(f'{HERE}/Fig{name}.png'); im.resize((im.size[0] // 3, im.size[1] // 3)).save(f'{HERE}/Fig{name}_prev.png')
def letters(fig, WMM, HMM, L):
    for s, x, y in L: fig.text(x / WMM, y / HMM, s, fontsize=9, fontweight='bold', va='top', ha='left')
GREY, RED = '#7F7F7F', OI['verm']

# ---------------- S1: how the drying time is counted
if 'S1' in WANT:
    X = pd.read_csv(D + 'figS1_example.csv', parse_dates=['date']); EVX = pd.read_csv(D + 'figS1_example_events.csv', parse_dates=['date'])
    CO = pd.read_csv(D + 'figS1_composite.csv'); HI = pd.read_csv(D + 'figS1_hist.csv', index_col=0)
    HI.loc[15] = HI.loc[15] + HI.loc[16]; HI = HI.drop(16)                          # same cap (15 days) as Figs 2b and 3d
    WMM, HMM = 183, 62
    fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
    a1 = ax_(12, 25, 78, 33); a2 = ax_(12, 8, 78, 14, sharex=a1); b = ax_(104, 8, 36, 50); c = ax_(152, 8, 30, 50)
    X = X[(X.date.dt.month >= 4) & (X.date.dt.month <= 10)]
    for _, e in EVX.iterrows():
        a1.axvspan(e.date + pd.Timedelta(hours=12), e.date + pd.Timedelta(days=int(e.rec), hours=12), color=OI['sky'], alpha=0.35, lw=0)
    a1.plot(X.date, X.ffmc, color='k', lw=0.8); a1.axhline(85, color='#4D4D4D', lw=0.6, ls='--')
    a1.set_ylim(20, 100); a1.set_yticks([40, 60, 85, 100]); a1.set_ylabel('FFMC', labelpad=1)
    a1.text(pd.Timestamp('1979-07-01'), 82.5, 'Dry fuels: FFMC ≥ 85', color='#4D4D4D', fontsize=6.5, va='top', ha='center')
    a1.text(pd.Timestamp('1979-07-15'), 32, 'Shading: drying time after each rain ≥ 10 mm', fontsize=6.5, ha='center', va='bottom', color='#3B7FA6')
    plt.setp(a1.get_xticklabels(), visible=False)
    a2.bar(X.date, X.rain, width=1.0, color=OI['blue']); a2.set_ylim(0, 45); a2.set_yticks([0, 20, 40]); a2.set_ylabel('Rain\n(mm)', labelpad=1)
    a2.set_xticks(pd.date_range('1979-04-01', '1979-10-01', freq='MS')); a2.set_xticklabels(['Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct'])
    a2.set_xlim(X.date.iloc[0], X.date.iloc[-1] + pd.Timedelta(days=1))
    for late, col, lb in [(0, GREY, '1979–2000'), (1, RED, '2001–2022')]:
        k = CO[CO.late == late]; b.fill_between(k.lag, k.lo, k.hi, color=col, alpha=0.25, lw=0); b.plot(k.lag, k.ffmc, color=col, lw=0.9, marker='o', ms=2)
        b.text(0.97, 0.08 + 0.09 * late, lb, transform=b.transAxes, color=col, fontsize=6.5, ha='right', va='bottom')
    b.axhline(85, color='k', lw=0.5, ls='--'); b.set_xlim(-3.5, 10.5); b.set_xticks([-3, 0, 3, 6, 9]); b.set_xticklabels(['−3', '0', '3', '6', '9'])
    b.set_ylim(50, 92); b.set_xlabel('Days from rain', labelpad=1); b.set_ylabel('Mean FFMC', labelpad=1)
    ib = b.inset_axes([0.47, 0.27, 0.49, 0.30])
    for late, col in [(0, GREY), (1, RED)]:
        k = CO[(CO.late == late) & (CO.lag >= 3)]; ib.plot(k.lag, k.ffmc, color=col, lw=0.8, marker='o', ms=1.5)
    ib.set_xlim(2.6, 10.4); ib.set_ylim(80, 86); ib.set_xticks([3, 6, 9]); ib.set_yticks([80, 85]); ib.tick_params(labelsize=6, length=1.5, pad=1)
    w = np.arange(len(HI)) + 1
    c.bar(w - 0.2, 100 * HI['0'], width=0.4, color=GREY); c.bar(w + 0.2, 100 * HI['1'], width=0.4, color=RED)
    c.set_xticks([1, 5, 10, 15]); c.set_xticklabels(['1', '5', '10', '≥15']); c.set_xlim(0.3, 15.7)
    c.set_xlabel('Drying time (days)', labelpad=1); c.set_ylabel('Share of rains (%)', labelpad=1)
    c.text(0.97, 0.97, '1979–2000', color=GREY, transform=c.transAxes, ha='right', va='top', fontsize=6.5)
    c.text(0.97, 0.88, '2001–2022', color=RED, transform=c.transAxes, ha='right', va='top', fontsize=6.5)
    letters(fig, WMM, HMM, [('a', 0.5, HMM - 0.5), ('b', 93, HMM - 0.5), ('c', 142, HMM - 0.5)])
    save(fig, 'S1')

# ---------------- S2: maps for alternative inputs and the Nelson model
if 'S2' in WANT:
    M = pd.read_csv(D + 'figS2_maps.csv')
    BOUNDS = [-30, -20, -15, -10, -5, -2, 2, 5, 10, 15, 20, 30]; CM, NORM = binned('BrBG', BOUNDS)
    P = [('era5', 'ERA5 weather and rain (reference)'), ('chirps', 'CHIRPS rain (1981–2000 vs 2001–2022)'), ('txrn', 'Daily maximum temperature, minimum humidity'),
         ('wind', 'Daily ERA5 noon wind'), ('nelson1', 'Nelson model, 1-h fuel sticks'), ('nelson10', 'Nelson model, 10-h fuel sticks')]
    WMM = 183; mw = 90; mh = mw * ASPECT; gap = 6.5; HMM = 3 * (mh + gap) + 12
    fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
    L = []
    for i, (k, t) in enumerate(P):
        r, cidx = divmod(i, 2); x = 1 + cidx * (mw + 1); y = HMM - (r + 1) * (mh + gap) + 1
        ax = ax_(x, y, mw, mh, projection=ROB); v = M[['lat', 'lon', k]].dropna()
        cellmap(ax, v.lat, v.lon, v[k].values, CM, NORM)
        ax.text(0.5, 1.02, t, transform=ax.transAxes, ha='center', va='bottom', fontsize=7)
        ax.text(0.06, 0.30, f'Median {m(v[k].median())}%\nDecline in {(v[k] < 0).mean() * 100:.0f}%\nof {len(v):,} cells', transform=ax.transAxes, fontsize=6.5, ha='left', va='top', path_effects=HALO)
        L.append(('abcdef'[i], x - 0.5, y + mh + gap - 0.5))
    hcolorbar(fig, ax_(46, 7.5, 91, 2.2), CM, NORM, BOUNDS, 'Change in drying time after rain, 2001–2022 vs 1979–2000 (%)')
    letters(fig, WMM, HMM, [(s, x, min(y, HMM - 0.5)) for s, x, y in L]); save(fig, 'S2')

# ---------------- S4: removing the temperature and humidity trends
if 'S4' in WANT:
    T = pd.read_csv(D + 'figS3_maps.csv')
    BOUNDS = [-30, -20, -15, -10, -5, -2, 2, 5, 10, 15, 20, 30]; CM, NORM = binned('BrBG', BOUNDS)
    TB = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2]; CT, NT = binned('YlOrRd', TB); CT.set_under(OI['sky'])   # cooling in blue
    WMM = 183; mw = 90; mh = mw * ASPECT; gap = 6.5; HMM = 2 * (mh + gap) + 23
    fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
    P = [('dT', 'Change in fire-season noon temperature'), ('warm', 'Humidity trend removed'), ('dry', 'Temperature trend removed'), ('neither', 'Both trends removed')]
    L = []
    for i, (k, t) in enumerate(P):
        r, cidx = divmod(i, 2); x = 1 + cidx * (mw + 1); y = HMM - (r + 1) * (mh + gap) + 1 - (10.5 if r == 1 else 0)
        ax = ax_(x, y, mw, mh, projection=ROB); v = T[['lat', 'lon', k]].dropna()
        cellmap(ax, v.lat, v.lon, v[k].values, *((CT, NT) if k == 'dT' else (CM, NORM)))
        ax.text(0.5, 1.02, t, transform=ax.transAxes, ha='center', va='bottom', fontsize=7)
        s = f'Median {m(v[k].median(), ".2f")} K\n{(v[k] > 0).mean() * 100:.0f}% warmer' if k == 'dT' else f'Median {m(v[k].median())}%\nDecline in {(v[k] < 0).mean() * 100:.0f}%\nof {len(v):,} cells'
        ax.text(0.06, 0.30, s, transform=ax.transAxes, fontsize=6.5, ha='left', va='top', path_effects=HALO)
        L.append(('abcd'[i], x - 0.5, y + mh + gap - 0.5))
        if i == 0: hcolorbar(fig, ax_(x + 8, y - 5.2, mw - 16, 2.2), CT, NT, TB, 'Temperature change (K)')
    hcolorbar(fig, ax_(46, 7.5, 91, 2.2), CM, NORM, BOUNDS, 'Change in drying time after rain, 2001–2022 vs 1979–2000 (%)')
    letters(fig, WMM, HMM, [(s, x, min(y, HMM - 0.5)) for s, x, y in L]); save(fig, 'S4')

# ---------------- S3: fire activity by FFMC class, dose-response slopes by region, and the visibility-corrected regional test
if 'S3' in WANT:
    import json, cartopy.crs as ccrs
    from figmaps import PC
    Q = pd.read_csv(D + 'figS4_ffmc_bins.csv'); G = pd.read_csv(D + 'figS3_slopes.csv'); G = G[G.n >= 300].sort_values('slope').reset_index(drop=True)
    FR = json.load(open(D + 'cloudfree_results.json')); PCL = pd.read_csv(D + 'figS3_cloudfree_cells.csv')
    WMM, HMM = 183, 198
    fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
    Y1 = 139                                                                   # bottom of the first row
    a = ax_(12, Y1 + 11, 52, 45); b = ax_(80, Y1 + 14.5, 102, 41.5)
    xs = np.arange(len(Q)); lab = [f'{int(l)}–{int(h)}' if h < 101 else f'≥{int(l)}' for l, h in zip(Q.lo, Q.hi)]; lab[0] = '<60'
    a.bar(xs, Q.ratio, color=[OI['sky'] if h <= 85 else OI['verm'] for h in Q.hi], width=0.75)
    a.axhline(1, color='k', lw=0.5, ls='--'); a.set_xticks(xs); a.set_xticklabels(lab, rotation=90, fontsize=6.5)
    a.set_ylim(0, 1.3); a.set_ylabel('Fire detections relative to\nthe same cell, year and month', labelpad=1)
    a.set_xlabel('FFMC', labelpad=1); a.tick_params(axis='x', pad=1)
    a.text(0.03, 0.97, 'Wet fuels', color='#2E8BC0', transform=a.transAxes, fontsize=6.5, va='top'); a.text(0.03, 0.88, 'Dry fuels', color=OI['verm'], transform=a.transAxes, fontsize=6.5, va='top')
    x = np.arange(len(G)); b.axhline(0, color='k', lw=0.5)
    b.bar(x, G.slope, color=['#BDBDBD' if s < 0 else '#6E6E6E' for s in G.slope], width=0.7)
    for i, r in G.iterrows(): b.plot([i, i], [r.lo, r.hi], color='k', lw=0.6, solid_capstyle='butt')
    b.set_xticks(x); b.set_xticklabels(G.ar6, rotation=90, fontsize=6.5); b.set_xlim(-0.7, len(G) - 0.3); b.tick_params(axis='x', pad=1)
    b.set_ylabel('Fewer fires per day of\ndrying time (days)', labelpad=1)
    b.text(0.02, 0.97, f'{(G.slope > 0).sum()} of {len(G)} regions positive, {(G.lo > 0).sum()} with intervals above zero', transform=b.transAxes, fontsize=6.5, va='top')
    # c: the five regions
    REGC = {'SAM': OI['green'], 'SES': OI['blue'], 'ESAF': OI['orange'], 'NAU': OI['purple'], 'CAF': OI['verm']}
    c = ax_(1, 100, 99, 32, projection=ccrs.Robinson(central_longitude=35)); c.set_extent([-84, 152, -40, 22], crs=PC); c.set_xticks([]); c.set_yticks([])
    import cartopy.feature as cfeature
    c.add_feature(cfeature.LAND.with_scale('110m'), facecolor='#E6E6E6', edgecolor='none', zorder=0)
    c.coastlines('110m', linewidth=0.3, color='#4D4D4D', zorder=2); c.spines['geo'].set_linewidth(0.5)
    for r_, g in PCL.groupby('ar6'):
        c.scatter(g.lon, g.lat, s=2.2, marker='s', color=REGC[r_], linewidths=0, transform=PC, zorder=1.5, rasterized=True)
    for r_, (lo, la, ha) in {'SAM': (-49, -2.5, 'right'), 'SES': (-46, -31, 'left'), 'CAF': (15, 17.5, 'center'), 'ESAF': (40, -24, 'left'),
                             'NAU': (133, -24, 'center')}.items():
        c.text(lo, la, r_, transform=PC, fontsize=6.5, color=REGC[r_], ha=ha, va='center', fontweight='bold', path_effects=HALO, zorder=3)
    # d, e: drying-time contrast and wet-fuel rate ratio for the five regions together, by period, in the checks and with each region left out
    get = lambda t: next(o for o in FR if o['tag'] == t)
    rows = [('five regions: main', 'Five regions, 2003–2022'), ('period 2003-2012', '2003–2012'), ('period 2013-2022', '2013–2022'),
            ('five regions: clear>=0.7', 'Clear sky ≥ 70%'), ('five regions: ten-day baseline', 'Ten-day baseline'), ('five regions: confidence 7-9', 'Fire confidence 7–9'),
            ('five regions: Terra', 'Terra only'), ('five regions: Aqua', 'Aqua only')] + [(f'without {r}', f'Without {r}') for r in ('SAM', 'SES', 'NAU', 'ESAF', 'CAF')]
    d = ax_(128, 60, 24, 72); e = ax_(158, 60, 24, 72)
    for ax, key, ci, ref in [(d, 'dT50', 'dT50_ci', 0), (e, 'M3_rr_wet', 'M3_ci', 1)]:
        ax.axvline(ref, color='k', lw=0.5, ls='--')
        for i_, (t, lb) in enumerate(rows):
            yy = len(rows) - 1 - i_; o = get(t); mc = 'k' if i_ == 0 else '#4D4D4D'
            lo_, hi_ = o[ci]
            if ax is d and hi_ > 9: hi_ = 8.8; ax.plot(8.8, yy, '>', ms=3, color=mc, mew=0)  # interval runs beyond the axis
            ax.plot([lo_, hi_], [yy, yy], color=mc, lw=0.7, solid_capstyle='butt'); ax.plot(o[key], yy, 'o', ms=3.4 if i_ == 0 else 3, color=mc)
        for yb in (len(rows) - 1.5, len(rows) - 3.5, len(rows) - 8.5): ax.axhline(yb, color='#BFBFBF', lw=0.5)
        ax.set_ylim(-0.6, len(rows) - 0.4); ax.set_yticks(range(len(rows) - 1, -1, -1)); ax.tick_params(axis='y', length=0)
        ax.set_yticklabels([lb for _, lb in rows] if ax is d else [], fontsize=6.5)
    d.set_xlim(-1, 9); d.set_xticks([0, 4, 8])
    e.set_xscale('log'); e.set_xlim(0.2, 1.25); e.set_xticks([0.2, 0.5, 1]); e.set_xticklabels(['0.2', '0.5', '1']); e.minorticks_off()
    d.set_xlabel('ΔT50, slow − fast (days)', labelpad=1); e.set_xlabel('Wet-fuel rate ratio', labelpad=1)
    # f: gradient-boosting check; g: five regions with clear-sky fractions of at least 70%; h-l: each region
    import json as _j
    SF = _j.load(open(D + 'fire_summary_modis.json')); CIF = _j.load(open(D + 'fig2_ci_modis.json'))
    f = ax_(13, 60, 26, 30)
    f.bar([0, 1], [SF['gbm_obs_over_pred_wet'], SF['gbm_obs_over_pred_dry']], color=[OI['sky'], OI['verm']], width=0.6)
    for xx, v in [(0, SF['gbm_obs_over_pred_wet']), (1, SF['gbm_obs_over_pred_dry'])]: f.text(xx, v + 0.03, m(v, '.2f'), ha='center', va='bottom', fontsize=6.5)
    f.axhline(1, color='k', lw=0.5, ls='--'); f.set_xlim(-0.6, 1.6); f.set_ylim(0, 1.45); f.set_xticks([0, 1]); f.set_xticklabels(['Wet\nfuels', 'Dry\nfuels'], fontsize=6.5)
    f.set_yticks([0, 0.5, 1]); f.set_ylabel('Observed / predicted\nfire detections', labelpad=1)
    f.text(0.03, 0.97, f'Wet/dry {m(CIF["gbm_ratio"], ".2f")}\n({m(CIF["gbm_ci"][0], ".2f")}–{m(CIF["gbm_ci"][1], ".2f")})', transform=f.transAxes, fontsize=6.5, va='top')
    CC = [OI['orange'], OI['sky'], OI['blue']]; T = np.arange(-10, 21)
    def sm(x): x = np.asarray(x, float); x = np.where(np.isfinite(x), x, np.nanmean(x)); return np.convolve(x, np.ones(3) / 3, 'same')
    def curves(ax, o, ttl, ylab, legend, col='k'):
        A = np.array(o['A'], float)
        ax.axhline(1, color='#4D4D4D', lw=0.5, ls='--'); ax.axhline(0.5, color='#4D4D4D', lw=0.5, ls=':'); ax.axvline(0, ymax=0.6, color='#4D4D4D', lw=0.5)
        jj = (T >= -7) & (T <= 15)
        for cl in range(3):
            ax.plot(T[jj], sm(A[cl])[jj], color=CC[cl], lw=1.0, label=['1 day', '2 days', '3–15 days'][cl])
            if o['T50'][cl] < 21: ax.plot(o['T50'][cl], 0.5, '|', ms=5, mew=1.2, color=CC[cl])
        ax.set_xlim(-7, 15); ax.set_xticks([-5, 0, 5, 10, 15]); ax.set_xticklabels(['−5', '0', '5', '10', '15'])
        ax.set_ylim(0, 1.95); ax.set_yticks([0, 0.5, 1, 1.5]); ax.set_xlabel('Days from rain', labelpad=1)
        if ylab: ax.set_ylabel('Fire activity relative to\npre-rain dry days', labelpad=1)
        else: ax.set_yticklabels([])
        lo, hi = o['dT50_ci']; cen = ' (slow class not at half by day 20)' if max(o['T50']) >= 21 else ''
        ax.text(0.03, 0.97, ttl, transform=ax.transAxes, fontsize=6.5, ha='left', va='top', color=col, fontweight='bold')
        ax.text(0.03, 0.86, f'ΔT50 {m(o["dT50"], ".1f")} ({m(lo, ".1f")} to {m(hi, ".1f")})', transform=ax.transAxes, fontsize=6.5, ha='left', va='top')
        if legend: ax.legend(loc='upper right', frameon=False, title='Drying time', title_fontsize=6.5, fontsize=6.5, handlelength=1.2, borderaxespad=0.1, labelspacing=0.2)
    g_ = ax_(54, 60, 58, 30); curves(g_, get('five regions: clear>=0.7'), 'Five regions, clear sky ≥ 70%', True, True)
    for k, r in enumerate(['SAM', 'SES', 'NAU', 'ESAF', 'CAF']):
        ax = ax_(13 + k * 34.2, 9, 31.5, 34); curves(ax, get(f'region {r}'), r, k == 0, False, REGC[r])
    letters(fig, WMM, HMM, [('a', 0.5, HMM - 0.5), ('b', 68, HMM - 0.5), ('c', 0.5, 134.5), ('d', 103, 134.5), ('e', 153, 134.5),
                            ('f', 0.5, 93), ('g', 42, 93), ('h', 0.5, 47.5), ('i', 45, 47.5), ('j', 79, 47.5), ('k', 113, 47.5), ('l', 147.5, 47.5)]); save(fig, 'S3')

# ---------------- S6: DAMIP forced changes per model and the loss without them
if 'S6' in WANT:
    RH = pd.read_csv(D + 'figS6_damip_TH.csv'); FX = pd.read_csv(D + 'fig4_damip.csv')        # v1.2: like-for-like monthly/daily-mean, fixed events
    obs = RH[RH.kind.isna()].set_index('source'); oT, oH = obs.loc['ERA5 daily mean', 'dT'], obs.loc['ERA5 daily mean', 'dRH']
    nT, nH = obs.loc['ERA5 noon', 'dT'], obs.loc['ERA5 noon', 'dRH']
    M = RH[RH.kind.notna()].pivot_table(index=['source', 'kind'], columns='var', values='val').reset_index().rename(columns={'source': 'model'})
    F5 = M.merge(FX[['model', 'kind', 'med', 'med_ref']], on=['model', 'kind'], how='outer'); F5['cf'] = F5.med / 100; F5['obs'] = F5.med_ref / 100
    mods = sorted(F5.model.unique()); short = {}
    WMM, HMM = 183, 66
    fig = plt.figure(figsize=(WMM * MM, HMM * MM)); ax_ = lambda x, y, w, h, **k: fig.add_axes([x / WMM, y / HMM, w / WMM, h / HMM], **k)
    a = ax_(11, 24, 38, 40); b = ax_(64, 24, 38, 40); c = ax_(118, 24, 64, 40)
    y = np.arange(len(mods))
    for ax, col, ov, lb in [(a, 'dT', oT, 'Fire-season warming (K)'), (b, 'dRH', oH, 'Fire-season humidity\nchange (%)')]:
        for j, (kind, cc) in enumerate([('ant', 'k'), ('ghg', OI['purple']), ('aer', OI['blue'])]):
            k = F5[F5.kind == kind].set_index('model').reindex(mods)[col]
            ax.scatter(np.arange(len(mods)) + (j - 1) * 0.22, k, s=12, color=cc, zorder=3, edgecolor='none')
        ax.axhline(ov, color=OI['verm'], lw=0.8); ax.axhline(nT if col == 'dT' else nH, color=OI['verm'], lw=0.6, ls='--'); ax.axhline(0, color='k', lw=0.4)
        ax.set_xticks(y); ax.set_xticklabels([short.get(q, q) for q in mods], rotation=90, fontsize=6.5); ax.set_xlim(-0.6, len(mods) - 0.4)
        ax.set_ylabel(lb, labelpad=1); ax.tick_params(axis='x', pad=1)
    a.set_ylim(-0.4, 1.2); b.set_ylim(-3.1, 0.6)
    a.text(0.03, 0.97, 'Observed: daily mean (solid), noon (dashed)', color=OI['verm'], transform=a.transAxes, fontsize=6, va='top')
    b.text(0.97, 0.03, 'Observed: daily mean,\nnoon (dashed)', color=OI['verm'], transform=b.transAxes, fontsize=6, ha='right', va='bottom')
    for t, cc, yy in [('Human, all', 'k', 0.64), ('Greenhouse gases', OI['purple'], 0.57), ('Aerosols', OI['blue'], 0.50)]:
        b.text(0.03, yy, t, color=cc, transform=b.transAxes, fontsize=6.5, va='top')
    kinds = [('ant', 'Human T and RH', 'k'), ('ant_T', 'Human, T only', '#6E6E6E'), ('ant_H', 'Human, RH only', '#BDBDBD'), ('ghg', 'Greenhouse gases', OI['purple']), ('aer', 'Aerosols', OI['blue'])]
    wb = 0.16
    for j, (kd, lb, cc) in enumerate(kinds):
        k = F5[F5.kind == kd].set_index('model').reindex(mods)
        c.bar(y + (j - 2) * wb, 100 * k.cf, width=wb, color=cc, label=lb)
    c.axhline(100 * F5.obs.median(), color=OI['verm'], lw=0.8); c.axhline(0, color='k', lw=0.4)
    c.set_xticks(y); c.set_xticklabels([short.get(q, q) for q in mods], rotation=90, fontsize=6.5); c.set_xlim(-0.6, len(mods) - 0.4); c.tick_params(axis='x', pad=1)
    c.set_ylim(-14, 0); c.set_yticks([-12, -8, -4, 0]); c.set_yticklabels(['−12', '−8', '−4', '0'])
    c.set_ylabel('Change in drying time with\nthe forced response removed (%)', labelpad=1)
    c.plot([], [], color=OI['verm'], lw=0.8, label='Observed')
    c.legend(loc='lower left', ncol=3, fontsize=6, frameon=False, handlelength=1, columnspacing=0.8, bbox_to_anchor=(0, -0.01))
    letters(fig, WMM, HMM, [('a', 0.5, HMM - 0.5), ('b', 51, HMM - 0.5), ('c', 104, HMM - 0.5)]); save(fig, 'S6')
