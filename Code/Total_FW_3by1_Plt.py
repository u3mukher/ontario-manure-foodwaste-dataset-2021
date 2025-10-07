"""
3 × 1 panels: Food-waste, FW-N, FW-P  (kg ha⁻¹ yr⁻¹)

Reads from:
  Ontario_30m_Manure_FW_N_P_kg_ha_fill-9999_OUTSIDE-ONTARIO.nc

Notes:
- Variables are already in kg/ha/yr.
- -9999 is used OUTSIDE Ontario; we mask it to NaN for plotting.
- ZEROS are masked to white (only non-zero values plotted), like your original FW script.
- For the FW MASS panel only, values are divided by 1000 to show ×10³ kg/ha/yr.
"""

# ------------------------------------------------------------------
# imports
# ------------------------------------------------------------------
import math
import xarray as xr
import rioxarray
import geopandas as gpd
import rasterio
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as mpatches
from rasterio.enums import Resampling

# ------------------------------------------------------------------
# paths & frame
# ------------------------------------------------------------------
NC_FILE   = "Ontario_30m_Dataset.nc"
BOUND_SHP = "Ontario_Boundary_5070.shp"

XMIN_3857, XMAX_3857 = -9_700_000, -8_000_000
YMIN_3857, YMAX_3857 =  5_050_000,  5_780_000

# nodata marker in the NC file (outside Ontario)
NODATA_VAL = -9999.0

# ------------------------------------------------------------------
# bin edges (unchanged)
# ------------------------------------------------------------------
BINS_FW   = [0, 0.6, 1.2, 3, 5, 12, 16]           # ×10³ kg/ha/yr
BINS_FW_N = [0, 5, 20, 50, 80, 120, 170]          # kg/ha/yr
BINS_FW_P = [0, 1, 4, 7, 10, 15, 23]              # kg/ha/yr

# ------------------------------------------------------------------
# interval labels  (<x,  x–y,  …)
# ------------------------------------------------------------------
_fmt = lambda v: f'{v:.1f}' if isinstance(v, float) and v % 1 else f'{int(v)}'

def interval_labels(edges):
    """Convert bin edges to labels like ['<0.6','0.6–1.2', …]."""
    first = f'<{_fmt(edges[1])}'
    spans = [f'{_fmt(lo)}–{_fmt(hi)}' for lo, hi in zip(edges[1:-1], edges[2:])]
    return [first] + spans

labels_fw   = interval_labels(BINS_FW)
labels_fw_n = interval_labels(BINS_FW_N)
labels_fw_p = interval_labels(BINS_FW_P)

# ------------------------------------------------------------------
# fonts – one master knob
# ------------------------------------------------------------------
FS_SCALE      = 1
TICK_FS       = int(11 * FS_SCALE)
LEG_TITLE_FS  = int(10 * FS_SCALE)
LETTER_FS     = int(14 * FS_SCALE)
SCALEBAR_FS   = int(10 * FS_SCALE)

# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------
def ensure_spatial_dims(da):
    """
    Make sure the DataArray has canonical x/y dims for rioxarray.
    Tries to map lon_5070->x and lat_5070->y (or infer from common names),
    writes CRS=EPSG:5070, and registers spatial dims.
    """
    dims = dict(da.sizes)
    xdim = None
    ydim = None
    for d in dims:
        dl = d.lower()
        if dl in ("lon_5070", "lon", "longitude", "x_5070", "x", "easting"):
            xdim = d if xdim is None else xdim
        if dl in ("lat_5070", "lat", "latitude", "y_5070", "y", "northing"):
            ydim = d if ydim is None else ydim
    if xdim is None or ydim is None:
        for d in dims:
            dl = d.lower()
            if xdim is None and ("lon" in dl or dl == "x" or "east" in dl):
                xdim = d
            if ydim is None and ("lat" in dl or dl == "y" or "north" in dl):
                ydim = d
    if xdim is None or ydim is None:
        raise ValueError(f"Could not infer spatial dims from {list(dims.keys())}")

    rename_map = {}
    if xdim != "x": rename_map[xdim] = "x"
    if ydim != "y": rename_map[ydim] = "y"
    if rename_map:
        da = da.rename(rename_map)

    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    da = da.rio.write_crs("EPSG:5070", inplace=False)
    return da

def downscale_self(da, max_pix=3000, nodata=None):
    """
    Downscale within current CRS using average resampling.
    Preserves nodata during reprojection.
    """
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
                      width=.10, height_frac=.80, xoff=-.10, gap_frac=0.05):
    """Vertical legend made of square patches with gaps and a two-line title."""
    pos = ax.get_position()
    h_tot = pos.height * height_frac
    cax = fig.add_axes([pos.x1 + xoff,
                        pos.y0 + (pos.height - h_tot)/2,
                        width, h_tot])
    cax.set_axis_off()

    n = len(colors)
    h_gap   = gap_frac
    h_patch = (1 - (n-1)*h_gap) / n

    for i, (c, lbl) in enumerate(zip(colors, labels)):
        y = 1 - i*(h_patch + h_gap) - h_patch
        cax.add_patch(
            mpatches.Rectangle((0.0, y), 0.4, h_patch,
                               facecolor=c, edgecolor='none',
                               transform=cax.transAxes)
        )
        cax.text(0.45, y + h_patch/2, lbl,
                 va='center', ha='left', fontsize=TICK_FS,
                 transform=cax.transAxes)

    fig.text(cax.get_position().x0 + width/2,
             cax.get_position().y1 + 0.01,
             '\n'.join(title_lines),
             ha='center', va='bottom', fontsize=LEG_TITLE_FS)

def add_scalebar_compass(ax, length_km=100,
                         x_pad_frac=0.3, y_pad_frac=0.3,
                         gap_multiplier=1.20):
    dx, dy = XMAX_3857 - XMIN_3857, YMAX_3857 - YMIN_3857
    length_m = length_km * 1000
    x_end   = XMAX_3857 - x_pad_frac*dx
    x_start = x_end - length_m
    y_bar   = YMIN_3857 + y_pad_frac*dy
    ax.plot([x_start, x_end], [y_bar, y_bar], color='k', linewidth=2)
    ax.text((x_start+x_end)/2, y_bar-0.02*dy, f'{length_km} km',
            ha='center', va='top', fontsize=SCALEBAR_FS)
    base = 0.0805*dy
    g    = 0.02*dy*gap_multiplier
    ax.annotate('', xy=(x_start+length_m/2, y_bar+g+base),
                xytext=(x_start+length_m/2, y_bar+g),
                arrowprops=dict(facecolor='black', edgecolor='black',
                                linewidth=2, headwidth=8, headlength=10))
    ax.text(x_start+length_m/2, y_bar+g+base+0.012*dy,
            'N', ha='center', va='bottom',
            fontsize=SCALEBAR_FS, fontweight='bold')

# ------------------------------------------------------------------
# data
# ------------------------------------------------------------------
bound = gpd.read_file(BOUND_SHP).to_crs(epsg=3857)
ds    = xr.open_dataset(
    NC_FILE,
    chunks={"lat_5070": 4096, "lon_5070": 4096},
    mask_and_scale=False
)

def prep_var(var, *, divide_by_1000_for_fw=False):
    """
    Prepare a variable from the NC:
      - Map to x/y, write CRS=5070
      - Write nodata = -9999 (do NOT convert to NaN yet)
      - Downscale in 5070 (average), preserving nodata
      - Reproject to 3857 (average), preserving nodata
      - Clip to Ontario boundary
      - Convert nodata to NaN (transparent in plots)
      - Mask zeros to NaN (white in plots), like original FW script
      - For FW MASS only, divide by 1000 to show ×10³ kg/ha/yr
    """
    da = ds[var].squeeze(drop=True).astype("float32")
    da = ensure_spatial_dims(da)

    # Declare nodata but don't yet convert to NaN
    da = da.rio.write_nodata(NODATA_VAL, inplace=False)

    # Downscale in native CRS (EPSG:5070)
    da = downscale_self(da, max_pix=3000, nodata=NODATA_VAL)

    # Reproject to Web Mercator
    da = da.rio.reproject("EPSG:3857", resampling=Resampling.average, nodata=NODATA_VAL)

    # Clip to Ontario AFTER reprojection
    da = da.rio.clip(bound.geometry, bound.crs, drop=False)

    # Convert nodata to NaN (transparent)
    da = da.where(da != NODATA_VAL)

    # Mask zeros to white (plot only non-zero values)
    da = da.where(da > 0)

    # For total FW mass only: convert kg/ha/yr to ×10³ kg/ha/yr
    if divide_by_1000_for_fw:
        da = da / 1000.0

    da.name = ""
    return da

# Variables likely named like:
#   "food_waste", "food_waste_N", "food_waste_P"
da_fw   = prep_var("food_waste",   divide_by_1000_for_fw=True)
da_fw_n = prep_var("food_waste_N")
da_fw_p = prep_var("food_waste_P")

# ------------------------------------------------------------------
# colour schemes (grey + 5 swatches)
# ------------------------------------------------------------------
grey = "#e0e0e0"

colors_fw   = [grey,  # Oranges (5-class)
               "#feedde", "#fdbe85", "#fd8d3c", "#e6550d", "#a63603"]

colors_fw_n = [grey,  # BuGn (5-class)
               "#edf8fb", "#b2e2e2", "#66c2a4", "#2ca25f", "#006d2c"]

colors_fw_p = [grey,  # YlOrBr (5-class)
               "#ffffd4", "#fee391", "#fec44f", "#fe9929", "#d95f0e"]

norm_fw   = BoundaryNorm(BINS_FW,   len(colors_fw),   clip=False)
norm_fw_n = BoundaryNorm(BINS_FW_N, len(colors_fw_n), clip=False)
norm_fw_p = BoundaryNorm(BINS_FW_P, len(colors_fw_p), clip=False)

layers = [
    (da_fw,   colors_fw,   norm_fw,   labels_fw,
     ["Total Food Waste Mass", "×10³ kg/ha/yr"]),
    (da_fw_n, colors_fw_n, norm_fw_n, labels_fw_n,
     ["Food Waste Nitrogen", "kg N/ha/yr"]),
    (da_fw_p, colors_fw_p, norm_fw_p, labels_fw_p,
     ["Food Waste Phosphorus", "kg P/ha/yr"])
]

# ------------------------------------------------------------------
# figure
# ------------------------------------------------------------------
ratio = (YMAX_3857 - YMIN_3857) / (XMAX_3857 - XMIN_3857)
fig_w = 6.5
fig_h = 3 * ratio * fig_w

fig, axes = plt.subplots(3, 1, figsize=(fig_w, fig_h), dpi=500)
plt.subplots_adjust(hspace=0.05)

# ------------------------------------------------------------------
# plot
# ------------------------------------------------------------------
label_x, label_y = 0.25, 0.50
letter = ord('a')

for idx, (ax, layer) in enumerate(zip(axes, layers)):
    da, colors, norm, labels, title_lines = layer
    cmap = ListedColormap(colors)
    cmap.set_bad(color=(1, 1, 1, 0))  # NaNs (nodata + zeros) transparent → white bg
    ax.set_facecolor('white')

    da.plot.imshow(ax=ax, cmap=cmap, norm=norm,
                   add_colorbar=False, add_labels=False,
                   interpolation="nearest")
    bound.boundary.plot(ax=ax, edgecolor='black', linewidth=0.35)

    ax.set_xlim(XMIN_3857, XMAX_3857)
    ax.set_ylim(YMIN_3857, YMAX_3857)
    ax.set_aspect('equal')
    ax.axis('off')

    # panel label
    ax.text(label_x, label_y, f"({chr(letter)})",
            transform=ax.transAxes,
            fontsize=LETTER_FS, fontweight='bold',
            va='center', ha='left')
    letter += 1

    # legend
    add_square_legend(ax, colors=colors, labels=labels,
                      title_lines=title_lines, fig=fig)

    # scale-bar & north arrow on top plot only
    if idx == 0:
        add_scalebar_compass(ax)

# ------------------------------------------------------------------
# save / show
# ------------------------------------------------------------------
fig.savefig("Total_FoodWaste_FW_N_P_3x1.png",
            dpi=1000, bbox_inches="tight", facecolor="white")
plt.show()





