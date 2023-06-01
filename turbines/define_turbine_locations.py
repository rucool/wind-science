#!/usr/bin/env python

"""
Author: Lori Garzio on 5/12/2023
Last modified: 6/1/2023
Define 2km grid spacing for simulated wind turbines in the NY/NJ wind lease areas
"""

import os
import xarray as xr
import numpy as np
import csv
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from shapely.geometry import Point
import functions.common as cf
import functions.plotting as pf
import cool_maps.plot as cplt
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def main(domain, save_file):
    if domain == '1km_ctrl':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_ctrl_processed/WRF_4.1_1km_Control_Processed_Dataset_Best'
    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=ds.time.values[0])
    #extent = [-74.6, -73.8, 38.85, 39.7]  # smaller region for debugging
    extent = [-75.2, -72.2, 38.1, 40.5]

    data_sub, lon, lat = cf.subset_wrf_grid(extent, ds.SST)

    # set up the map
    kwargs = dict()
    kwargs['oceancolor'] = 'none'
    kwargs['coast'] = 'high'
    fig, ax = cplt.create(extent, **kwargs)

    lease = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/BOEM_Renewable_Energy_Shapefiles_2_2023/Wind_Lease_Outlines_2_2023.shp'
    bwargs = dict()
    bwargs['edgecolor'] = 'magenta'
    bwargs['zorder'] = 10
    bwargs['linewidth'] = .6
    geoms = pf.map_add_boem_outlines(ax, lease, **bwargs)

    # # plot all points
    # ax.scatter(lon, lat, s=.1, transform=ccrs.PlateCarree())
    #
    # plt.savefig(os.path.join(os.path.dirname(save_file), 'wrf_1km_gridpoints.png'), dpi=200)
    plt.close()

    mask = np.ones(np.shape(lon), dtype=bool)
    lon_array = np.array([], dtype='float32')
    lat_array = np.array([], dtype='float32')

    for ii, geom in enumerate(geoms):
        for i, loni in enumerate(lon):
            for j, lonj in enumerate(lon[i]):
                if geom.contains(Point(lon[i, j], lat[i, j])):  # find grid points within WEAs
                    if np.logical_and(i % 2 == 0, j % 2 == 0):  # find every other grid point/line
                        mask[i, j] = False
                        lon_array = np.append(lon_array, lon[i, j].values)
                        lat_array = np.append(lat_array, lat[i, j].values)

    lon.values[mask] = np.nan
    lat.values[mask] = np.nan

    # plot points within the lease areas
    fig, ax = cplt.create(extent, **kwargs)
    geoms = pf.map_add_boem_outlines(ax, lease, **bwargs)
    ax.scatter(lon, lat, s=.1, transform=ccrs.PlateCarree())

    plt.savefig(os.path.join(os.path.dirname(save_file), 'wrf_1km_gridpoints_leaseareas_2km.png'), dpi=200)
    plt.close()

    # export tab-delimited .txt file of turbine locations
    turb_type = np.ones(np.shape(lon_array), dtype=int)
    txt_file = os.path.join(os.path.dirname(save_file), 'windturbines.txt.2kmspacing.allleaseareas')
    with open(txt_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['lat', 'lon', 'turb_type'], delimiter=' ')
        for a, b, c in zip(lat_array, lon_array, turb_type):
            writer.writerow({'lat': a, 'lon': b, 'turb_type': c})


if __name__ == '__main__':
    domain = '1km_ctrl'
    savefile = '/Users/garzio/Documents/rucool/bpu/wrf/windturbs/plots/wrf_1km_turbines.png'
    main(domain, savefile)
