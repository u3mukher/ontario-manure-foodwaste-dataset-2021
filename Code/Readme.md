### `memory_efficient_query.ipynb`

This Jupyter notebook shows how to open the Ontario manure & food-waste NetCDF efficiently and extract just the data you need, without loading the full 30 m grid into memory. In particular, it will:

1. Open the NetCDF with Dask-backed Xarray (`chunks="auto"`).  
2. Select a single variable (e.g. `total_manure_N`).  
3. Subset that variable over a user-specified spatial window.  
4. Compute the percentage of cells with non-zero values (i.e. spatial coverage).  
5. Plot the subset to visualize spatial patterns.

### Running the example notebook

1. Download `memory_efficient_query.ipynb` into the `Codes/` folder.
2. Download `Ontario_30m_Manure_FW_N_P_kg_ha.nc` into the `Data/` folder.
3. Open the notebook in Jupyter (e.g. `jupyter lab Codes/memory_efficient_query.ipynb`) and run the cells to see how to query a single variable efficiently.

