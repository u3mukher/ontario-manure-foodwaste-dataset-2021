### `memory_efficient_query.ipynb`

This Jupyter notebook shows how to open the Ontario manure & food waste NetCDF in a **memory-efficient** way using chunked reads via `xarray`. Rather than loading the entire 30 m grid at once, it reads only the data needed for distribution statistics. The notebook demonstrates:

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

