#!/usr/bin/env python

"""
Author: Lori Garzio on 9/25/2024
Last modified: Lori Garzio 9/25/2024
Take the subset files generated from wrf_data_wrangler_grid.py and generate averages for the user-specified time range
"""

import os
import xarray as xr
import pandas as pd
import datetime as dt
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console


def main(start_str, end_str, ver, v):
    # combine files into one dataset
    files = [
        f'/Users/garzio/Documents/rucool/Miles/RMI/2024_Ocean_Mixing/data/{ver}/ruwrf_{ver}_{v}_20220601_20220630.nc',
        f'/Users/garzio/Documents/rucool/Miles/RMI/2024_Ocean_Mixing/data/{ver}/ruwrf_{ver}_{v}_20220701_20220731.nc',
        f'/Users/garzio/Documents/rucool/Miles/RMI/2024_Ocean_Mixing/data/{ver}/ruwrf_{ver}_{v}_20220801_20220831.nc']

    for i, f in enumerate(files):
        if i == 0:
            ds = xr.open_dataset(f)
        else:
            ds = xr.concat([ds, xr.open_dataset(f)], dim="time")

    # select time range
    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)
    ds = ds.sel(time=slice(start_date, end_date))

    if v in ['ws', 'ws10']:
        # calculate wind speed
        if v == 'ws':
            uvector = ds['U']
            vvector = ds['V']
        else:
            uvector = ds['U10']
            vvector = ds['V10']

        ds[v] = cf.wind_uv_to_spd(uvector, vvector)
        ds[v].attrs['units'] = uvector.attrs['units']
        ds[v].attrs['long_name'] = 'Wind Speed'

    varname = f'{v}_avg'
    ds[varname] = ds[v].mean(dim='time')
    ds[varname].attrs['units'] = ds[v].attrs['units']
    ds[varname].attrs['comment'] = f'{v} average from {start_str} to {end_str}'
    try:
        ds[varname].attrs['long_name'] = f'Average {ds[v].attrs["long_name"]}'
    except KeyError:
        ds[varname].attrs['long_name'] = f'Average {v}'

    # Add compression to all variables
    encoding = {}
    for k in ds.data_vars:
        encoding[k] = {'zlib': True, 'complevel': 1}

    save_file = os.path.join(os.path.dirname(f), f'ruwrf_{ver}_{v}_avg_{start_str}_{end_str}.nc')
    ds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4')


if __name__ == '__main__':
    start = '20220601'
    end = '20220831'
    version = '1km_wf2km_nyb'  # 1km_ctrl 1km_wf2km_nyb
    variable = 'ws10'  # HFX Q2 TKE_PBL UST T2 TEMP ws10 ws
    main(start, end, version, variable)
