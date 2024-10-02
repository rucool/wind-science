#!/usr/bin/env python

"""
Author: Lori Garzio on 9/30/2024
Last modified: 10/1/2024
Compare SST for NOW23 and RU-WRF datasets with user-specified NDBC buoy
"""

import numpy as np
import pandas as pd
import glob
import os
import datetime as dt
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def find_closest_model_point(da, mlon, mlat, blon, blat):
    # calculate the sum of the absolute value distance between the model location and buoy location
    a = abs(mlat - blat) + abs(mlon - blon)

    # find the indices of the minimum value in the array calculated above
    i, j = np.unravel_index(a.argmin(), a.shape)

    return(da[:, i, j])


def format_date_axis(axis):
    #datef = mdates.DateFormatter('%m-%d\n%H:%M')
    datef = mdates.DateFormatter('%b-%d\n%Y')
    axis.xaxis.set_major_formatter(datef)


def get_buoy_sst(buoy_name, all_years, t0, t1):
    # get buoy SST for entire time range
    buoy_dict = {'t': np.array([], dtype='datetime64[ns]'), 'sst': np.array([]), 'sst_units': '',
                 'lon': '', 'lat': ''}
    for y in all_years:
        buoy_fname = f'{buoy_name}h{str(y)}'
        buoydap = f'https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/{buoy_name}/{buoy_fname}.nc'
        try:
            ds = xr.open_dataset(buoydap)
            ds = ds.sel(time=slice(t0, t1))
            if len(ds['time']) > 0:
                ds_sst = np.squeeze(ds['sea_surface_temperature'].values)
                ds_sst[ds_sst == ds.sea_surface_temperature.encoding['_FillValue']] = np.nan  # convert fill values to nans
                buoy_dict['t'] = np.append(buoy_dict['t'], ds['time'].values)
                buoy_dict['sst'] = np.append(buoy_dict['sst'], ds_sst)
                buoy_dict['sst_units'] = ds['sea_surface_temperature'].units
                buoy_dict['lon'] = ds['longitude'].values[0]
                buoy_dict['lat'] = ds['latitude'].values[0]
        except OSError:
            print('File does not exist: {}'.format(buoydap))

    return buoy_dict


def main(start_str, end_str, filedir, save_dir, buoy, v):
    os.makedirs(save_dir, exist_ok=True)

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    dr = pd.date_range(start_date, end_date)
    years = np.unique(dr.year)

    # grab buoy data
    buoy_data = get_buoy_sst(buoy, years, start_date, end_date)

    # grab now23 dataset
    if v == 'sst':
        now23var = 'surface_sea_temperature'
        wrfvar = 'SST'
        ylab = 'SST (\N{DEGREE SIGN}C)'
        title = f'SST comparison at NDBC Buoy {buoy}'
    elif v == 'skin':
        now23var = 'skin_temperature'
        wrfvar = 'TSK'
        ylab = 'Skin Temperature (\N{DEGREE SIGN}C)'
        title = f'Skin Temperature comparison at NDBC Buoy {buoy}'

    now23_dict = {'t': np.array([], dtype='datetime64[ns]'), v: np.array([])}
    files = sorted(glob.glob(os.path.join(filedir, '*.nc')))
    for f in files:
        ds = xr.open_dataset(f)
        ds = ds.sel(time=slice(start_date, end_date))
        sst = find_closest_model_point(ds[now23var], ds.lon, ds.lat, buoy_data['lon'], buoy_data['lat'])
        now23_dict['t'] = np.append(now23_dict['t'], sst['time'].values)
        now23_dict[v] = np.append(now23_dict[v], sst)

    # get ru-wrf data
    ruwrf_dict = {'t': np.array([], dtype='datetime64[ns]'), v: np.array([])}
    ds = xr.open_dataset('https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best')
    ds = ds.sel(time=slice(start_date, end_date))
    sst = find_closest_model_point(ds[wrfvar], ds.XLONG, ds.XLAT, buoy_data['lon'], buoy_data['lat'])
    ruwrf_dict['t'] = np.append(ruwrf_dict['t'], sst['time'].values)
    ruwrf_dict[v] = np.append(ruwrf_dict[v], sst - 273.15)

    fig, ax = plt.subplots(figsize=(13, 7))

    ax.plot(buoy_data['t'], buoy_data['sst'], color='k', label=f'Buoy {buoy}')
    ax.plot(now23_dict['t'], now23_dict[v], color='#d95f02', label='NOW23')  # orange
    ax.plot(ruwrf_dict['t'], ruwrf_dict[v], color='#7570b3', label='RU-WRF 3km')  # purple

    ax.set_ylim([18, 28])
    ax.set_ylabel(ylab)
    ax.set_xlabel('Time (GMT)')
    ax.legend(loc='best', fontsize=10)
    format_date_axis(ax)

    plt.title(title)

    save_file = f'{v}_now23_vs_wrf_vs_buoy{buoy}_{start_str}-{end_str}.png'
    plt.savefig(os.path.join(save_dir, save_file), dpi=200)
    plt.close()


if __name__ == '__main__':
    start = '20200801'
    end = '20200831'
    now23files = '/Users/garzio/Documents/rucool/bpu/now23_comparison/now23_data/august_2020_test'
    save_dir = '/Users/garzio/Documents/rucool/bpu/now23_comparison/sst_comparison/buoy_comparison_timeseries'
    buoy = '44091'
    variable = 'skin'  # sst skin
    main(start, end, now23files, save_dir, buoy, variable)
