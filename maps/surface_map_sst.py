#!/usr/bin/env python

"""
Author: Lori Garzio on 3/14/2023
Last modified: 3/14/2023
Creates a surface map of sea surface temperature for a user-defined date, domain, and region/extent
"""

import argparse
import sys
import numpy as np
import pandas as pd
import glob
import os
import datetime as dt
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import cmocean as cmo
import cool_maps.plot as cplt
import functions.plotting as pf
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def subset_grid(ext, data, lon_name, lat_name):
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


def main(args):
    ymd = args.ymd
    domain = args.domain
    plot_region = args.plot_region
    clims = args.clims
    save_dir = args.save_dir

    ymd_dt = dt.datetime.strptime(ymd, '%Y%m%d')
    ym = ymd[0:6]

    save_dir = os.path.join(save_dir, 'surface_maps_sst', ym)
    os.makedirs(save_dir, exist_ok=True)

    plt_region = pf.plot_regions()
    extent = plt_region[plot_region]['extent']

    if domain == '3km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
        title_label = '3km'
    elif domain == '1km_wf2km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_wf2km_processed/WRF_4.1_1km_with_Wind_Farm_Processed_Dataset_Best'
        title_label = '1km Wind Farm'
    elif domain == '1km_ctrl':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_ctrl_processed/WRF_4.1_1km_Control_Processed_Dataset_Best'
        title_label = '1km Control'
    else:
        raise ValueError('Invalid domain specified')

    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=ymd_dt + dt.timedelta(hours=1))

    sst = ds.SST - 273.15  # convert from K to degrees C

    sst_sub, lon, lat = subset_grid(extent, sst, 'XLONG', 'XLAT')

    kwargs = dict()
    kwargs['oceancolor'] = 'none'
    fig, ax = cplt.create(extent, **kwargs)
    contour_list = [5, 10, 15, 20, 25, 30]
    pf.add_contours(ax, lon, lat, sst_sub.values, contour_list)

    # define color map
    bins = clims[1] - clims[0]
    cmap = cmo.cm.thermal
    levels = MaxNLocator(nbins=bins).tick_values(clims[0], clims[1])
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    kwargs = dict()
    kwargs['norm_clevs'] = norm
    kwargs['extend'] = 'both'
    kwargs['cmap'] = cmap
    kwargs['clab'] = 'SST (\N{DEGREE SIGN}C)'
    pf.plot_pcolormesh(fig, ax, lon, lat, sst_sub.values, **kwargs)
    ax.set_title(f'{title_label} SST {pd.to_datetime(ymd).strftime("%Y-%m-%d")}', pad=8)

    #lease = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/BOEM-Renewable-Energy-Shapefiles_11_2_2022/Wind_Lease_Outlines_11_2_2022.shp'
    lease = glob.glob('/home/coolgroup/bpu/mapdata/shapefiles/BOEM-Renewable-Energy-Shapefiles-current/Wind_Lease_Outlines*.shp')[0]
    kwargs = dict()
    kwargs['edgecolor'] = 'magenta'
    pf.map_add_boem_outlines(ax, lease, **kwargs)

    save_file = f'ruwrf_sst_{domain}_{plot_region}_{ymd}.png'
    plt.savefig(os.path.join(save_dir, save_file), dpi=200)
    plt.close()


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description='Plot WRF SST',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('ymd',
                            type=str,
                            help='Year-month-day to plot in the format YYYYmmdd (e.g. 20220101.')

    arg_parser.add_argument('-d', '--domain',
                            dest='domain',
                            default='1km_ctrl',
                            type=str,
                            choices=['3km', '1km_wf2km', '1km_ctrl'],
                            help='Operational: 3km, research 1km with simulated windfarm: 1km_wf2km,'
                                 'research 1km control: 1km_ctrl')

    arg_parser.add_argument('-plot_region',
                            default='windturb',
                            type=str,
                            choices=['full_grid', 'mab', 'nj', 'southern_nj', 'windturb'],
                            help='Region to plot')

    arg_parser.add_argument('-clims',
                            default=[18, 31],
                            type=list,
                            help='Colorbar limits [min, max]')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/windturbs/plots',
                            type=str,
                            help='Full directory path to save output plots.')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
