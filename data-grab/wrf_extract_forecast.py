#!/usr/bin/env python

"""
Author: Lori Garzio on 10/12/2023
Last modified: Lori Garzio 10/12/2023
Grab RU-WRF forecast data for a user-specified domain, date and location
"""

import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console


def main(fd, dom, loc):
    start_date = dt.datetime.strptime(fd, '%Y%m%d')

    if dom == '3km':
        mlink = 'http://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_fmrc.ncd'
    elif dom == '9km':
        mlink = 'http://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_9km_processed/WRF_4.1_9km_Processed_Dataset_fmrc.ncd'
    else:
        raise ValueError('Invalid domain specified')

    ds = xr.open_dataset(mlink)

    # select data for the specified run date
    ds = ds.sel(run=start_date)

    # select data at the specified location
    # calculate the sum of the absolute value distance between the model location and buoy location
    lat = ds.XLAT
    lon = ds.XLONG
    a = abs(lat - loc[0]) + abs(lon - loc[1])

    # find the indices of the minimum value in the array calculated above
    # subset the coordinate dimensions using those indices
    i, j = np.unravel_index(a.argmin(), a.shape)
    ds = ds.isel(south_north=i, west_east=j)


if __name__ == '__main__':
    forecast_date = '20181201'
    domain = '9km'  # 3km 9km
    location = [42.65, -73.77]  # Albany NY [lat, lon]
    main(forecast_date, domain, location)
