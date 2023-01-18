#!/usr/bin/env python

"""
Author: Lori Garzio on 11/12/2022
Last modified: 11/12/2022
Generate a map of the lease areas and their centers
"""

import pandas as pd
import matplotlib.pyplot as plt
import functions.hurricanes_plotting as hp
import cartopy.crs as ccrs
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def main(save_file, extent):
    # set up the map
    kwargs = dict()
    kwargs['zoom_coastline'] = False
    fig, ax = hp.map_create(extent, **kwargs)

    lease = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/BOEM-Renewable-Energy-Shapefiles_11_2_2022/Wind_Lease_Outlines_11_2_2022.shp'
    plan = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/BOEM-Renewable-Energy-Shapefiles_11_2_2022/Wind_Planning_Area_Outlines_11_2_2022.shp'
    kwargs = dict()
    kwargs['edgecolor'] = 'dimgray'
    hp.map_add_boem_outlines(ax, lease, **kwargs)

    # kwargs['edgecolor'] = 'lightgray'
    # hp.map_add_boem_outlines(ax, plan, **kwargs)

    lease_cent_csv = '/Users/garzio/Documents/repo/rucool/wind-science/files/lease_centroids.csv'
    loc_df = pd.read_csv(lease_cent_csv)
    for i, row in loc_df.iterrows():
        if row.state == 'New Jersey':
            ax.scatter(row['long'], row['lat'], s=20, marker='o', c='cyan', edgecolor='k', transform=ccrs.PlateCarree(),
                       zorder=20)
        if row.state == 'NY/NJ':
            ax.scatter(row['long'], row['lat'], s=20, marker='o', c='magenta', edgecolor='k', transform=ccrs.PlateCarree(),
                       zorder=20)

    plt.savefig(save_file, dpi=200)
    plt.close()


if __name__ == '__main__':
    savefile = '/Users/garzio/Documents/rucool/bpu/wrf/capacity_factor_analysis/map_lease_areas.png'
    #extent = [-74.9, -72.3, 38.7, 40.3]  # zoomed in to WEA
    extent = [-75.5, -72, 38.4, 40.6]
    main(savefile, extent)
