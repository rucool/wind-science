#!/usr/bin/env python

"""
Author: Lori Garzio on 6/28/2022
Last modified: 1/18/2023
Plot surface map of the WRF validation locations
"""

import pandas as pd
import ast
import functions.plotting as pf
import matplotlib.pyplot as plt
import cool_maps.plot as cplt
import cartopy.crs as ccrs
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def add_validation_points(axis, text, longitude, latitude, toffset):
    axis.scatter(longitude, latitude, s=40, marker='o', c='magenta', edgecolor='k', transform=ccrs.PlateCarree(), zorder=20)

    lon_offset = toffset[0]
    lat_offset = toffset[1]

    text_lon = longitude + lon_offset
    text_lat = latitude + lat_offset
    axis.text(text_lon, text_lat, text, fontsize=10, transform=ccrs.PlateCarree(), zorder=20)


def main(save_file, shape_files, ndbc, zoom):
    if zoom:
        extent = [-74.9, -72.3, 38.7, 40.3]  # zoomed in to WEA
    else:
        extent = [-75.5, -72, 38.4, 40.6]

    # set up the map
    fig, ax = cplt.create(extent)

    if shape_files:
        #lease = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/BOEM-Renewable-Energy-Shapefiles_11_2_2022/Wind_Lease_Outlines_11_2_2022.shp'
        #plan = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/BOEM-Renewable-Energy-Shapefiles_11_2_2022/Wind_Planning_Area_Outlines_11_2_2022.shp'  #find -name Wind_Lease_Outlines_11_2_2022.shp
        lease = '/users/jameskim/Documents/rucool/bpu/maps/BOEM_Renewable_Energy_Shapefiles_1/Wind_Lease_Outlines_2_2023.shp'
        plan = '/users/jameskim/Documents/rucool/bpu/maps/BOEM_Renewable_Energy_Shapefiles_1/BOEM_Wind_Planning_Areas_04_2023.shp'


        kwargs = dict()
        kwargs['edgecolor'] = 'dimgray'
        pf.map_add_boem_outlines(ax, lease, **kwargs)

        kwargs['edgecolor'] = 'lightgray'
        pf.map_add_boem_outlines(ax, plan, **kwargs)

    df = pd.read_csv('/users/jameskim/Documents/rucool/Repositories/wind-science/files/wrf_validation_points.csv')

    for i, row in df.iterrows():
        if row['name'] == 'RUOYC':
            continue
        if zoom:
            offset = ast.literal_eval(row['text_offset_zoom'])
        else:
            offset = ast.literal_eval(row['text_offset'])
        add_validation_points(ax, row['name'], row['longitude'], row['latitude'], offset)

    if ndbc:
        df = pd.read_csv('/users/jameskim/Documents/rucool/Repositories/wind-science/files/ndbc_buoys.csv')
        for i, row in df.iterrows():
            if zoom:
                offset = ast.literal_eval(row['text_offset_zoom'])
            else:
                offset = ast.literal_eval(row['text_offset'])
            add_validation_points(ax, row['name'], row['longitude'], row['latitude'], offset)

    plt.savefig(save_file, dpi=200)
    plt.close()


if __name__ == '__main__':
    savefile = '/Users/jameskim/Documents/rucool/bpu/maps/wrf_validation_points_20230118.png'
    shpfiles = True  # True False
    add_ndbc = False  # True False
    zoom_wea = True
    main(savefile, shpfiles, add_ndbc, zoom_wea)
