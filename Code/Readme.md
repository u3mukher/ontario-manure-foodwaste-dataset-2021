### `memory_efficient_query.ipynb`

This Jupyter notebook shows how to open the Ontario manure & food waste NetCDF in a memory-efficient way using chunked reads via `xarray`. Rather than loading the entire 30 m grid at once, it reads only the data needed for distribution statistics. The notebook demonstrates:

1. Opening the NetCDF file with `xarray` using chunked access (`chunks="auto"`).  
2. Extracting a single variable of interest (e.g., `total_manure`).  
3. Masking out zero-value cells on the fly to focus analysis on meaningful data.  
4. Computing key distribution metrics across all non-zero cells:  
   - 25th percentile (Q1)  
   - Median (50th percentile)  
   - 75th percentile (Q3)  
   - Minimum non-zero value  
   - Maximum value  
5. Performing all calculations without ever loading the full grid into memory.

This example is intended for users new to NetCDF or `xarray` who need to derive summary statistics from large spatial datasets on standard hardware.

---
### `manure_plot.py` — 3×1 Map Panels (Manure, N, P)

This script generates a three-panel (3×1) map visualization showing total manure mass and associated nutrient contents (N and P) across Ontario.  
It is designed for large rasters and uses **xarray** + **rioxarray** for efficient reprojection and nodata handling.

#### Demonstrates
- Opening a large NetCDF dataset with `xarray` and maintaining geospatial metadata using `rioxarray`.  
- Mapping coordinate names to standard x/y, assigning `CRS = EPSG:5070`, and reprojecting to `EPSG:3857`.  
- Downscaling and reprojection with nodata-aware averaging to reduce memory and preserve edge accuracy.  
- Clipping to the Ontario boundary and masking the `–9999` nodata cells to transparent for clean borders.  
- Rendering three vertically stacked panels:
  1. **Total Manure Mass** (×10³ kg ha⁻¹ yr⁻¹; lowest non-zero bin shown in grey)  
  2. **Manure Nitrogen** (kg N ha⁻¹ yr⁻¹)  
  3. **Manure Phosphorus** (kg P ha⁻¹ yr⁻¹)  
- Adding square-patch legends, a scale bar, and a north arrow directly on the plots.

#### Key Behaviors
- **Bins:** Controlled by `BINS_MAN`, `BINS_N`, and `BINS_P` arrays.  
- **Zero Handling:** Zeros are masked to white while labels keep their `<1`, `<5`, `<2` formatting.  
- **Nodata:** Outside-Ontario pixels (`–9999`) remain transparent after reprojection.  
- **Performance:** Large rasters are automatically downscaled (average resampling) before plotting.

#### Inputs
- `Ontario_30m_Dataset.nc` (or the kg ha⁻¹ version)  
- `Ontario_Boundary_5070.shp` (converted to EPSG 3857 inside the script)

#### Output
- `Total_Manure_N_P_3x1.png` — high-resolution figure suitable for publication.

#### Dependencies
`xarray`, `rioxarray`, `geopandas`, `rasterio`, `matplotlib`


---

### `food_waste_plot.py` — 3×1 Map Panels (FW, FW-N, FW-P)

This script renders a three-panel figure showing food-waste mass and nutrient intensities (N, P) per hectare, following the same structure as the manure map.

#### Demonstrates
- Reading food-waste variables from the NetCDF file using chunked `xarray` reads.  
- Normalizing coordinate names, writing `CRS = EPSG:5070`, downscaling, and reprojecting to `EPSG:3857`.  
- Clipping to Ontario boundaries and masking `–9999` (outside area) and zero-value pixels to white.  
- Creating stacked panels with categorical bins:
  1. **Total Food Waste Mass** (×10³ kg ha⁻¹ yr⁻¹ — divided by 1000 for readability)  
  2. **Food Waste Nitrogen** (kg N ha⁻¹ yr⁻¹)  
  3. **Food Waste Phosphorus** (kg P ha⁻¹ yr⁻¹)  
- Adding compact legends, a north arrow, and a 100 km scale bar.

#### Key Behaviors
- **Bins:** `BINS_FW`, `BINS_FW_N`, and `BINS_FW_P` define color intervals.  
- **Units:** Only the FW-mass panel is rescaled to ×10³ kg ha⁻¹ yr⁻¹.  
- **Zeros:** Masked to white; legend labels remain in `<0.6`, `<5`, `<1` form.  
- **Nodata:** `–9999` handled consistently through reprojection and clipping.

#### Inputs
- `Ontario_30m_Dataset.nc` (or the kg ha⁻¹ variant)  
- `Ontario_Boundary_5070.shp`

#### Output
- `Total_FoodWaste_FW_N_P_3x1.png` — high-resolution PNG output.

#### Dependencies
`xarray`, `rioxarray`, `geopandas`, `rasterio`, `matplotlib`
