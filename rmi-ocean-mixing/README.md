# wind-science rmi-ocean-mixing
Tools for grabbing a subset of RU-WRF data from [THREDDS](https://tds.marine.rutgers.edu/thredds/catalog/cool/ruwrf/catalog.html), calculating averages and generating surface maps to compare to the analysis outlined in Golbazi et al 2022 [DOI 10.1088/1748-9326/ac6e49](https://iopscience.iop.org/article/10.1088/1748-9326/ac6e49/meta).

See the installation instructions on the main page for cloning this repo and generating the wind-science environment.

## Installation Instructions
Add the channel conda-forge to your .condarc. You can find out more about conda-forge from their website: https://conda-forge.org/

`conda config --add channels conda-forge`

Clone the wind-science repository.

`git clone https://github.com/rucool/wind-science.git`

Change your current working directory to the location that you downloaded wind-science/rmi-ocean-mixing. 

`cd /Users/lgarzio/Documents/repo/wind-science/rmi-ocean-mixing`

Create conda environment from the included environment.yml file:

`conda env create -f environment.yml`

Once the environment is done building, activate the environment:

`conda activate rmi-ocean-mixing`

Install the toolbox to the conda environment from the root directory of the wind-science/rmi-ocean-mixing toolbox:

`pip install .`

The toolbox should now be installed to your conda environment.


## Scripts

1. [wrf_data_wrangler_grid.py](https://github.com/rucool/wind-science/blob/master/rmi-ocean-mixing/wrf_data_wrangler_grid.py): Grabs RU-WRF output for a user-defined model domain/version and time range for a subset of variables and heights, and saves each file to a user-defined directory. For this analysis, we grabbed summer 2022 data monthly (June, July, August) from model domain/version 1km_ctrl (1km domain with no turbines) and 1km\_wf2km\_nyb (1km domain with a simulated windfarm including turbines spaced 2km apart and filling all wind energy lease areas in the New York Bight as of early 2023). This results in 3 NetCDF files (2022-06-01 to 2022-06-30, 2022-07-01 to 2022-07-31, and 2022-08-01 to 2022-08-31) for each variable. It's necessary to grab the data in small chunks because of the volume of model output.
	- variables: T2, TEMP, U, V, U10, V10, UST, Q2, TKE_PBL, HFX, TSK
	- heights: surface, 120m, 160m, 200m, 300m
	
2. [data_averages.py](https://github.com/rucool/wind-science/blob/master/rmi-ocean-mixing/data_averages.py): Generate averages of the datasets subset using wrf\_data\_wrangler_grid.py for a user-specified time range. For this analysis, we calculated averages for the entire summer 2022 (2022-06-01 to 2022-08-31).

The resulting data files from steps 1-2 are available [here](https://marine.rutgers.edu/~lgarzio/rmi_ocean_mixing/) and can be used to generate the surface maps using the scripts below:

3. [surface\_maps_avgdiff.py](https://github.com/rucool/wind-science/blob/master/rmi-ocean-mixing/surface_maps_avgdiff.py): Generate pcolormesh surface maps for each variable at multiple heights (if applicable) from the files averaged using data_averages.py. Maps include:
	- Control model output
	- Model output with simulated turbines
	- Difference (turbines minus control)
	- Proportion reduction of turbines compared to the control (difference divided by control)

4. [surface\_map_tutorial.ipynb](https://github.com/rucool/wind-science/blob/master/rmi-ocean-mixing/surface_map_tutorial.ipynb): Jupyter notebook tutorial for generating surface maps from [files](https://marine.rutgers.edu/~lgarzio/rmi_ocean_mixing/) containing summer 2022 average model output. To run a Jupyter notebook you can either use a Python IDE that supports Jupyter (e.g. [VS Code](https://code.visualstudio.com/)), or in your terminal navigate to the location that you downloaded wind-science/rmi-ocean-mixing:

`cd /Users/lgarzio/Documents/repo/wind-science/rmi-ocean-mixing`

Activate your environment:

`conda activate rmi-ocean-mixing`

Type jupyter notebook:

`jupyter notebook`

This will start a Jupyter Server in your browser that serves those notebooks in your local directory. You can click on the surface_map_tutorial.ipynb in your browser window to run through the tutorial to generate your own maps. To shut down the Jupyter server, close the windows in your browser and type Ctrl+C in your terminal.

## Acknowledgement
This work was supported by the New Jersey Research and Monitoring Initiative (RMI): New Jersey Board of Public Utilities (NJBPU) and New Jersey Department of Environmental Protection (NJDEP)
