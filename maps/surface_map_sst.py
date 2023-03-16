#!/usr/bin/env python

"""
Author: Lori Garzio on 3/14/2023
Last modified: 3/16/2023
Creates a surface map of sea surface temperature for a user-defined date, domain, and region/extent. If a date range
is specified, average SST is plotted.
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
    start_str = args.start
    end_str = args.end
    domain = args.domain
    plot_region = args.plot_region
    clims = args.clims
    save_dir = args.save_dir

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d')

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

    if end_date - start_date == dt.timedelta(0):
        save_dir = os.path.join(save_dir, 'surface_maps_sst', 'daily', start_str[0:6])
        save_file = f'ruwrf_sst_{domain}_{plot_region}_{start_str}.png'
        title = f'{title_label} SST {start_date.strftime("%Y-%m-%d")}'
        date_slice = start_date + dt.timedelta(hours=1)
        color_label = 'SST (\N{DEGREE SIGN}C)'
    else:
        save_dir = os.path.join(save_dir, 'surface_maps_sst', 'averages', start_str[0:4])
        save_file = f'ruwrf_avg_sst_{domain}_{plot_region}_{start_str}_{end_str}.png'
        title = f'{title_label} Average SST {start_date.strftime("%Y-%m-%d")} to {end_date.strftime("%Y-%m-%d")}'
        date_slice = pd.date_range(start_date, end_date) + dt.timedelta(hours=1)
        color_label = 'Average SST (\N{DEGREE SIGN}C)'

    os.makedirs(save_dir, exist_ok=True)

    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=date_slice)

    sst = ds.SST - 273.15  # convert from K to degrees C

    try:
        sst = sst.mean(dim='time')
    except ValueError:
        print('no averaging required')

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
    kwargs['clab'] = color_label
    pf.plot_pcolormesh(fig, ax, lon, lat, sst_sub.values, **kwargs)
    ax.set_title(title, pad=8)

    lease = glob.glob('/home/coolgroup/bpu/mapdata/shapefiles/BOEM-Renewable-Energy-Shapefiles-current/Wind_Lease_Outlines*.shp')[0]
    kwargs = dict()
    kwargs['edgecolor'] = 'magenta'
    pf.map_add_boem_outlines(ax, lease, **kwargs)

    plt.savefig(os.path.join(save_dir, save_file), dpi=200)
    plt.close()


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description='Plot WRF SST',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-s', '--start',
                            dest='start',
                            default='20220801',
                            type=str,
                            help='Start Date in format YYYYMMDD. If end date equals start date, a surface map for just'
                                 'that day is provided. Otherwise a surface map of average SST for the date range'
                                 'is provided.')

    arg_parser.add_argument('-e', '--end',
                            dest='end',
                            default='20220801',
                            type=str,
                            help='End Date in format YYYYMMDD. If end date equals start date, a surface map for just'
                                 'that day is provided. Otherwise a surface map of average SST for the date range'
                                 'is provided.')

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
