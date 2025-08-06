### `memory_efficient_query.ipynb`

This Jupyter notebook shows how to open the Ontario manure & food waste NetCDF in a **memory-efficient** way using chunked reads via `xarray`. Rather than loading the entire 30 m grid at once, it reads only the data needed for basic statistics. The notebook demonstrates:

1. Opening the NetCDF file with `xarray` using chunked access (`chunks="auto"`).  
2. Extracting a single variable of interest (e.g., `total_manure_N`).  
3. Computing the **Ontario-wide mean** and **non-zero spatial coverage**.  
4. Computing the same statistics for a specific region (Windsor → Ottawa extent).  
5. Avoiding full-grid reads while obtaining useful summary metrics.

This example is intended for users new to NetCDF or `xarray` who want to work with large spatial datasets without high memory use.
