# wind-science rmi-ocean-mixing
Tools for grabbing a subset of RU-WRF data from [THREDDS](https://tds.marine.rutgers.edu/thredds/catalog/cool/ruwrf/catalog.html), calculating averages and generating surface maps to compare to the analysis outlined in Golbazi et al 2022 [DOI 10.1088/1748-9326/ac6e49](https://iopscience.iop.org/article/10.1088/1748-9326/ac6e49/meta).

See the installation instructions on the main page for cloning this repo and generating the wind-science environment.

## Scripts

1. [wrf_data_wrangler_grid.py](https://github.com/rucool/wind-science/blob/master/rmi-ocean-mixing/wrf_data_wrangler_grid.py): Grabs RU-WRF output for a user-defined model domain/version and time range for a subset of variables and heights, and saves each file to a user-defined directory. For this analysis, we grabbed summer 2022 data monthly (June, July, August) from model domain/version 1km_ctrl (1km domain with no turbines) and 1km\_wf2km\_nyb (1km domain with a simulated windfarm including turbines spaced 2km apart and filling all wind energy lease areas in the New York Bight as of early 2023). This results in 3 NetCDF files (2022-06-01 to 2022-06-30, 2022-07-01 to 2022-07-31, and 2022-08-01 to 2022-08-31) for each variable. It's necessary to grab the data in small chunks because of the volume of model output.
	- variables: T2, TEMP, U, V, U10, V10, UST, Q2, TKE_PBL, HFX, TSK
	- heights: surface, 120m, 160m, 200m, 300m
	
2. [data_averages.py](https://github.com/rucool/wind-science/blob/master/rmi-ocean-mixing/data_averages.py): Generate averages of the datasets subset using wrf\_data\_wrangler_grid.py for a user-specified time range. For this analysis, we calculated averages for the entire summer 2022 (2022-06-01 to 2022-08-31).

3. [surface\_maps_avgdiff.py](https://github.com/rucool/wind-science/blob/master/rmi-ocean-mixing/surface_maps_avgdiff.py): Generate pcolormesh surface maps for each variable at multiple heights (if applicable) from the files averaged using data_averages.py. Maps include:
	- Control model output
	- Model output with simulated turbines
	- Difference (turbines minus control)
	- Proportion reduction of turbines compared to the control (difference divided by control)

## Acknowledgement
This work was supported by the New Jersey Research and Monitoring Initiative (RMI): New Jersey Board of Public Utilities (NJBPU) and New Jersey Department of Environmental Protection (NJDEP)
