#!/usr/bin/env python

"""
Author: Lori Garzio on 9/30/2024
Last modified: 10/1/2024
Creates a 2-panel surface map of sea surface temperature from RU-WRF 3km and now23 dataset
"""

import numpy as np
import pandas as pd
import glob
import os
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cmocean as cmo
import cool_maps.plot as cplt
import functions.plotting as pf
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def subset_grid(ext, data, lon_name, lat_name):
    if len(np.shape(data[lon_name])) == 1:
        lonx, laty = np.meshgrid(data[lon_name], data[lat_name])
    else:
        lonx = data[lon_name]
        laty = data[lat_name]

    lon_idx = np.logical_and(lonx > ext[0], lonx < ext[1])
    lat_idx = np.logical_and(laty > ext[2], laty < ext[3])

    # find i and j indices of lon/lat in boundaries
    ind = np.where(np.logical_and(lat_idx, lon_idx))

    # subset data from min i,j corner to max i,j corner
    # there will be some points outside of defined boundaries because grid is not rectangular
    data_sub = np.squeeze(data)[range(np.min(ind[0]), np.max(ind[0]) + 1), range(np.min(ind[1]), np.max(ind[1]) + 1)]
    lon = data_sub[lon_name]
    lat = data_sub[lat_name]

    return data_sub, lon, lat


def main(filedir, save_dir, plot_region, cmin, cmax, v):
    save_dir = os.path.join(save_dir, v)
    os.makedirs(save_dir, exist_ok=True)

    if plot_region == 'custom':
        extent = [-76.5, -68.5, 37, 42.2]
    else:
        plt_region = pf.plot_regions()
        extent = plt_region[plot_region]['extent']

    ds = xr.open_dataset('https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best')

    files = sorted(glob.glob(os.path.join(filedir, '*.nc')))
    for f in files:
        nds = xr.open_dataset(f)
        t0 = np.nanmin(nds.time.values)
        tf = np.nanmax(nds.time.values)

        for t in nds.time.values:
            nds1 = nds.sel(time=t)
            ds1 = ds.sel(time=t)

            tm = pd.to_datetime(t)
            tmsavestr = tm.strftime("%Y%m%dT%H%M")
            tmstr = tm.strftime("%Y-%m-%dT%H:%M")
            save_file = f'{v}_now23_vs_wrf_{tmsavestr}.png'
            if v == 'sst':
                main_title = f'NOW23 vs RU-WRF SST: {tmstr}'
                sst1 = nds1.surface_sea_temperature
                sst2 = ds1.SST
            elif v == 'skin':
                main_title = f'NOW23 vs RU-WRF Skin Temperature: {tmstr}'
                sst1 = nds1.skin_temperature
                sst2 = ds1.TSK

            sst1_sub, now23lon, now23lat = subset_grid(extent, sst1, 'lon', 'lat')
            sst2_sub, lon, lat = subset_grid(extent, sst2, 'XLONG', 'XLAT')

            # turn nans to zero
            sst1_sub = np.nan_to_num(sst1_sub)
            now23lon = np.nan_to_num(now23lon)
            now23lat = np.nan_to_num(now23lat)

            sst2_sub = sst2_sub - 273.15  # convert RU-WRF data from K to degrees C

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), sharey=True,
                                           subplot_kw=dict(projection=ccrs.Mercator()))
            fig.suptitle(main_title, fontsize=17, y=.93)

            kwargs = dict()
            kwargs['oceancolor'] = 'none'
            kwargs['landcolor'] = 'none'
            kwargs['decimal_degrees'] = True
            kwargs['coast'] = 'high'  # low full high
            kwargs['ax'] = ax1
            cplt.create(extent, **kwargs)
            kwargs['tick_label_left'] = False
            kwargs['ax'] = ax2
            cplt.create(extent, **kwargs)

            contour_list = [5, 10, 15, 20, 25, 30]
            pf.add_contours(ax1, now23lon, now23lat, sst1_sub, contour_list)
            pf.add_contours(ax2, lon, lat, sst2_sub.values, contour_list)

            # define color map
            bins = cmax - cmin
            cmap = cmo.cm.thermal
            levels = MaxNLocator(nbins=bins).tick_values(cmin, cmax)
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

            if v == 'sst':
                clab = 'SST (\N{DEGREE SIGN}C)'
                title1 = 'NOW23 SST'
                title2 = 'RU-WRF SST'
            elif v == 'skin':
                clab = 'Skin Temperature (\N{DEGREE SIGN}C)'
                title1 = 'NOW23 skin temp'
                title2 = 'RU-WRF TSK'

            kwargs = dict()
            kwargs['norm_clevs'] = norm
            kwargs['extend'] = 'both'
            kwargs['cmap'] = cmap
            kwargs['clab'] = clab
            pf.plot_pcolormesh(fig, ax1, now23lon, now23lat, sst1_sub, **kwargs)
            ax1.set_title(title1, pad=8)

            pf.plot_pcolormesh(fig, ax2, lon, lat, sst2_sub.values, **kwargs)
            ax2.set_title(title2, pad=8)

            plt.savefig(os.path.join(save_dir, save_file), dpi=200)
            plt.close()


if __name__ == '__main__':
    now23files = '/Users/garzio/Documents/rucool/bpu/now23_comparison/now23_data/august_2020_test'
    save_dir = '/Users/garzio/Documents/rucool/bpu/now23_comparison/sst_comparison/surface_maps'
    plot_region = 'custom'  # full_grid nj mab custom
    cmin = 18
    cmax = 29
    variable = 'skin'  # sst skin
    main(now23files, save_dir, plot_region, cmin, cmax, variable)
