"""
One 3×1 panel figure using:
  Ontario_30m_Manure_FW_N_P_kg_ha_fill-9999_OUTSIDE-ONTARIO.nc

Behavior:
- Data are already kg/ha/yr.
- -9999 outside-Ontario fill values are masked to transparent.
- ZEROS are masked to white (not plotted).
- First bin (grey) remains labeled as "<1", "<5", "<2" etc. (labels unchanged).
- Output: Total_Manure_N_P_3x1_nonzero_GREYfirst.png
"""

# ------------------------------------------------------------------
# imports
# ------------------------------------------------------------------
import math
import xarray as xr
import rioxarray
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as mpatches
from rasterio.enums import Resampling

# ------------------------------------------------------------------
# files & constants
# ------------------------------------------------------------------
NC_FILE   = "Ontario_30m_Dataset.nc"
BOUND_SHP = "Ontario_Boundary_5070.shp"

XMIN_3857, XMAX_3857 = -9_700_000, -8_000_000
YMIN_3857, YMAX_3857 =  5_050_000,  5_780_000

# ── bins (kg/ha/yr) ────────────────────────────────────────────────
BINS_MAN = [v * 1000 for v in [0, 1, 10, 20, 30, 40, 50]]
BINS_N   = [0, 5, 20, 60, 100, 180, 380]
BINS_P   = [0, 2, 7, 20, 40, 55, 67]

# ── fonts – one master knob ────────────────────────────────────────
FS_SCALE      = 1
TICK_FS       = int(11 * FS_SCALE)
LEG_TITLE_FS  = int(10 * FS_SCALE)
LETTER_FS     = int(14 * FS_SCALE)
SCALEBAR_FS   = int(10  * FS_SCALE)

# ── misc constants ────────────────────────────────────────────────
MAX_PIX     = 3000
BOUNDARY_LW = 0.35
NODATA_VAL  = -9999.0

VAR_MAN = "total_manure"
VAR_N   = "total_manure_N"
VAR_P   = "total_manure_P"

# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------
def ensure_spatial_dims(da):
    """
    Ensure x/y dims for rioxarray. Maps lon_5070->x and lat_5070->y
    (or infers from common names), writes CRS=EPSG:5070, and registers dims.
    """
    dims = dict(da.sizes)
    xdim = ydim = None
    for d in dims:
        dl = d.lower()
        if dl in ("lon_5070","lon","longitude","x_5070","x","easting") and xdim is None:
            xdim = d
        if dl in ("lat_5070","lat","latitude","y_5070","y","northing") and ydim is None:
            ydim = d
    if xdim is None or ydim is None:
        for d in dims:
            dl = d.lower()
            if xdim is None and ("lon" in dl or dl == "x" or "east" in dl):
                xdim = d
            if ydim is None and ("lat" in dl or dl == "y" or "north" in dl):
                ydim = d
    if xdim is None or ydim is None:
        raise ValueError(f"Could not infer spatial dims from {list(dims.keys())}")

    if xdim != "x" or ydim != "y":
        da = da.rename({xdim: "x", ydim: "y"})
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    da = da.rio.write_crs("EPSG:5070", inplace=False)
    return da

def downscale_self(da, max_pix=3000, nodata=None):
    """Downscale in-place CRS using average resampling (nodata-aware)."""
    ny, nx = da.sizes["y"], da.sizes["x"]
    f = max(1, math.ceil(max(ny, nx) / max_pix))
    if f == 1:
        return da.astype("float32")
    rx, ry = da.rio.resolution()
    return da.rio.reproject(
        da.rio.crs,
        resolution=(abs(rx) * f, abs(ry) * f),
        resampling=Resampling.average,
        nodata=nodata
    ).astype("float32")

def add_square_legend(ax, *, colors, labels, title_lines, fig,
                      width=.10, height_frac=.80, xoff=-.1,
                      gap_frac=0.05, edgecolors=None):
    """Vertical legend made of square patches with gaps and a two-line title."""
    pos = ax.get_position()
    h_total = pos.height * height_frac
    cax = fig.add_axes([pos.x1 + xoff,
                        pos.y0 + (pos.height - h_total) / 2,
                        width, h_total])
    cax.set_axis_off()

    n = len(colors)
    h_gap = gap_frac
    h_patch = (1 - (n - 1) * h_gap) / n

    for i, (color, lbl) in enumerate(zip(colors, labels)):
        edge = 'none' if edgecolors is None else edgecolors[i]
        y_top = 1 - i * (h_patch + h_gap) - h_patch
        cax.add_patch(
            mpatches.Rectangle((0.0, y_top), 0.4, h_patch,
                               facecolor=color, edgecolor=edge,
                               linewidth=1.0,
                               transform=cax.transAxes)
        )
        cax.text(0.45, y_top + h_patch / 2, lbl,
                 va='center', ha='left', fontsize=TICK_FS,
                 transform=cax.transAxes)

    fig.text(cax.get_position().x0 + width / 2,
             cax.get_position().y1 + 0.01,
             '\n'.join(title_lines),
             ha='center', va='bottom', fontsize=LEG_TITLE_FS)

def add_scalebar_compass(ax, length_km=100,
                         x_pad_frac=0.3, y_pad_frac=0.3,
                         gap_multiplier=1.20):
    dx, dy = XMAX_3857 - XMIN_3857, YMAX_3857 - YMIN_3857
    length_m = length_km * 1000
    x_end   = XMAX_3857 - x_pad_frac * dx
    x_start = x_end - length_m
    y_bar   = YMIN_3857 + y_pad_frac * dy

    ax.plot([x_start, x_end], [y_bar, y_bar], color='k', linewidth=2)
    ax.text((x_start + x_end)/2, y_bar - 0.02*dy,
            f'{length_km} km', ha='center', va='top', fontsize=SCALEBAR_FS)

    base_to_head = 0.0805*dy
    gap          = 0.02*dy*gap_multiplier
    ax.annotate('', xy=(x_start+length_m/2, y_bar+gap+base_to_head),
                xytext=(x_start+length_m/2, y_bar+gap),
                arrowprops=dict(facecolor='black', edgecolor='black',
                                linewidth=2, headwidth=8, headlength=10))
    ax.text(x_start+length_m/2, y_bar+gap+base_to_head+0.012*dy,
            'N', ha='center', va='bottom',
            fontsize=SCALEBAR_FS, fontweight='bold')

# ------------------------------------------------------------------
# data prep
# ------------------------------------------------------------------
bound = gpd.read_file(BOUND_SHP).to_crs(epsg=3857)
ds    = xr.open_dataset(
    NC_FILE,
    chunks={"lat_5070":4096, "lon_5070":4096},
    mask_and_scale=False
)

def prep_var(var):
    """
    Prepare variable:
      - set x/y dims, write nodata=-9999
      - downscale in EPSG:5070 (average), nodata-aware
      - reproject to EPSG:3857 (average), nodata-aware
      - clip to Ontario
      - mask nodata to NaN (transparent)
    """
    da = ds[var].squeeze(drop=True).astype("float32")
    da = ensure_spatial_dims(da)
    da = da.rio.write_nodata(NODATA_VAL, inplace=False)
    da = downscale_self(da, MAX_PIX, nodata=NODATA_VAL)
    da = da.rio.reproject("EPSG:3857", resampling=Resampling.average, nodata=NODATA_VAL)
    da = da.rio.clip(bound.geometry, bound.crs, drop=False)
    da = da.where(da != NODATA_VAL)
    da.name = ""
    return da

da_man = prep_var(VAR_MAN)
da_N   = prep_var(VAR_N)
da_P   = prep_var(VAR_P)

# ------------------------------------------------------------------
# colour schemes and norms
# ------------------------------------------------------------------
grey = "#e0e0e0"   # lighter universal grey

colors_man = [grey,
              "#ffcccc", "#ff9999", "#ff6666", "#ff3333", "#cc0000"]

colors_n   = [grey,
              "#e0f3db", "#a8ddb5", "#7bccc4", "#43a2ca", "#0868ac"]

colors_p   = [grey,
              "#ffffd4", "#fee391", "#fec44f", "#fe9929", "#d95f0e"]

norm_man = BoundaryNorm(BINS_MAN, len(colors_man), clip=False)
norm_n   = BoundaryNorm(BINS_N,   len(colors_n),   clip=False)
norm_p   = BoundaryNorm(BINS_P,   len(colors_p),   clip=False)

# Labels (keep original "<x" style even though zeros are masked)
labels_man = ['<1',  '1-10',  '10-20',  '20-30',  '30-40',  '40-50']
labels_n   = ['<5',  '5-20',  '20-60',  '60-100', '100-180', '180-380']
labels_p   = ['<2',  '2-7',   '7-20',   '20-40',  '40-55',  '55-67']

# ------------------------------------------------------------------
# ONLY PLOT 2: zeros masked to white; first bin is GREY for NON-ZERO values
# ------------------------------------------------------------------
# Mask zeros ONLY for plotting in this figure (nodata already NaN)
da_man_nz = da_man.where(da_man > 0)
da_N_nz   = da_N.where(da_N > 0)
da_P_nz   = da_P.where(da_P > 0)

# Figure layout
ratio = (YMAX_3857 - YMIN_3857) / (XMAX_3857 - XMIN_3857)
fig_w = 6.5
fig_h = 3 * ratio * fig_w

fig, axes = plt.subplots(3, 1, figsize=(fig_w, fig_h), dpi=500)
plt.subplots_adjust(hspace=0.05)
fig.patch.set_facecolor('white')

layers = [
    (da_man_nz, colors_man, norm_man, labels_man,
     ["Total Manure Mass", "×10³ kg/ha/yr"]),
    (da_N_nz,   colors_n,   norm_n,   labels_n,
     ["Manure Nitrogen", "kg N/ha/yr"]),
    (da_P_nz,   colors_p,   norm_p,   labels_p,
     ["Manure Phosphorus", "kg P/ha/yr"])
]

label_x, label_y = 0.25, 0.50
letter = ord('a')

for idx, (ax, layer) in enumerate(zip(axes, layers)):
    da, colors, norm, labels, title_lines = layer
    cmap = ListedColormap(colors)         # first color = grey (lowest non-zero bin)
    cmap.set_bad(color=(1, 1, 1, 0))      # NaNs (zeros & outside) transparent → white bg
    ax.set_facecolor('white')

    da.plot.imshow(ax=ax, cmap=cmap, norm=norm,
                   add_colorbar=False, add_labels=False,
                   interpolation="nearest")
    bound.boundary.plot(ax=ax, edgecolor='black', linewidth=BOUNDARY_LW)

    ax.set_xlim(XMIN_3857, XMAX_3857)
    ax.set_ylim(YMIN_3857, YMAX_3857)
    ax.set_aspect('equal'); ax.axis('off')

    ax.text(label_x, label_y, f"({chr(letter)})",
            transform=ax.transAxes, fontsize=LETTER_FS,
            fontweight='bold', va='center', ha='left')
    letter += 1

    add_square_legend(ax,
                      colors=colors,
                      labels=labels,
                      edgecolors=None,
                      title_lines=title_lines,
                      fig=fig)

    if idx == 0:
        add_scalebar_compass(ax)

fig.savefig("Total_Manure_N_P_3x1.png",
            dpi=1000, bbox_inches="tight", facecolor="white")

plt.show()
