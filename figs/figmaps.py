# Shared Robinson map helper for the 1.5-degree cell maps (same projection, extent and styling as the main figures).
import numpy as np, matplotlib
import cartopy.crs as ccrs, cartopy.feature as cfeature
PC = ccrs.PlateCarree(); ROB = ccrs.Robinson(central_longitude=10)
yproj = lambda lat: ROB.transform_points(PC, np.zeros_like(np.asarray(lat, float)), np.asarray(lat, float))[:, 1]
Y0, Y1 = yproj([-57, 80]); X0, X1 = ROB.x_limits
ASPECT = (Y1 - Y0) / (X1 - X0)
LON_E = np.arange(241) * 1.5 - 0.75; LAT_E = np.clip(np.arange(122) * 1.5 - 90.75, -90, 90)

def cellmap(ax, lat, lon, val, cmap, norm):
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_xticks([]); ax.set_yticks([])
    ax.add_feature(cfeature.LAND.with_scale('110m'), facecolor='#D9D9D9', edgecolor='none', zorder=0)
    g = np.full((121, 240), np.nan)
    ia = np.rint((np.asarray(lat) + 90) / 1.5).astype(int); io = np.rint((np.asarray(lon) % 360) / 1.5).astype(int) % 240
    g[ia, io] = val
    ax.pcolormesh(LON_E, LAT_E, g, cmap=cmap, norm=norm, transform=PC, zorder=1, rasterized=True)
    ax.coastlines('110m', linewidth=0.3, color='#4D4D4D', zorder=2); ax.spines['geo'].set_linewidth(0.5)
    ax.gridlines(crs=PC, linewidth=0.3, color='#B3B3B3', linestyle='--', xlocs=range(-180, 181, 60), ylocs=[-30, 0, 30, 60], zorder=1.5)

def hcolorbar(fig, cax, cmap, norm, bounds, label, fmt=lambda b: format(b, 'g').replace('-', '−')):
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal', extend='both', ticks=bounds)
    cb.ax.set_xticklabels([fmt(b) for b in bounds]); cb.ax.tick_params(labelsize=6.5, length=2, pad=1); cb.outline.set_linewidth(0.4)
    cb.set_label(label, fontsize=7, labelpad=1.5); return cb

def stipple(ax, lat, lon, mask, s=0.25):
    """Small dots on the cells where mask is true (e.g. at least 80% of models agree on the sign)."""
    lat = np.asarray(lat)[np.asarray(mask)]; lon = np.asarray(lon)[np.asarray(mask)]
    ax.scatter(lon, lat, s=s, color='#262626', linewidths=0, transform=PC, zorder=1.8, rasterized=True)
