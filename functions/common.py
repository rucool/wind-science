#! /usr/bin/env python

import logging
import numpy as np
import pandas as pd
import datetime as dt


def daterange_interval(interval, time_array):
    """
    Break up a time xarray data array into the plotting interval specified
    Options: monthly, seasonal
    """
    dr_interval = []
    if interval == 'monthly':
        groups = time_array.resample(time='MS', keep_attrs=True).groups
    elif interval == 'seasonal':
        groups = time_array.resample(time='QS-DEC', keep_attrs=True).groups
    else:
        raise ValueError(f'Invalid interval specified: {interval}')

    for group, values in groups.items():
        min_date = np.nanmin(time_array[values])
        max_date = np.nanmax(time_array[values])
        dr_interval.append([min_date, max_date])

    return dr_interval


def return_seabreeze_datetimes(csvfile='/home/wrfadmin/toolboxes/wind-science/files/BPU_Seabreeze.csv'):
    """
    Returns two datetime indices of seabreeze vs non-seabreeze days at hourly intervals (to match WRF output)
    csvfile: csv file containing seabreeze and non-seabreeze days, default is ../files/BPU_Seabreeze.csv file
    """
    df = pd.read_csv(csvfile)

    # get seabreeze datetimes
    df_sb = df[df['sea_breeze'] == 'Y']
    sb_dates = np.array(pd.to_datetime(df_sb['date']))
    sb_datetimes = [pd.date_range(pd.to_datetime(x), pd.to_datetime(x) + dt.timedelta(hours=23), freq='H') for x in
                    sb_dates]
    sb_datetimes = pd.to_datetime(sorted([inner for outer in sb_datetimes for inner in outer]))

    # get non-seabreeze datetimes
    df_nosb = df[df['sea_breeze'] == 'N']
    nosb_dates = np.array(pd.to_datetime(df_nosb['date']))
    nosb_datetimes = [pd.date_range(pd.to_datetime(x), pd.to_datetime(x) + dt.timedelta(hours=23), freq='H') for x
                      in nosb_dates]
    nosb_datetimes = pd.to_datetime(sorted([inner for outer in nosb_datetimes for inner in outer]))

    return sb_datetimes, nosb_datetimes


def return_upwelling_datetimes(csvfile='/home/wrfadmin/toolboxes/wind-science/files/NJUpwellingforBPU.csv'):
    """
    Returns two datetime indices of upwelling vs non-upwelling days at hourly intervals (to match WRF output)
    csvfile: csv file containing upwelling and non-upwelling days, default is ../files/NJUpwellingforBPU.csv file
    """
    df = pd.read_csv(csvfile)

    # drop May
    df = df[df.Month != 'May']

    # get seabreeze datetimes
    df['Month'] = df['Month'].map(lambda x: dt.datetime.strptime(x, '%B').month)
    df['date'] = pd.to_datetime(dict(year=df.Year, month=df.Month, day=df.Day))
    df_up = df[df['AVHRR'] == 1]
    up_dates = np.array(pd.to_datetime(df_up['date']))
    up_datetimes = [pd.date_range(pd.to_datetime(x), pd.to_datetime(x) + dt.timedelta(hours=23), freq='H') for x in
                    up_dates]
    up_datetimes = pd.to_datetime(sorted([inner for outer in up_datetimes for inner in outer]))

    # get non-seabreeze datetimes
    df_noup = df[df['AVHRR'] == 0]
    noup_dates = np.array(pd.to_datetime(df_noup['date']))
    noup_datetimes = [pd.date_range(pd.to_datetime(x), pd.to_datetime(x) + dt.timedelta(hours=23), freq='H') for x
                      in noup_dates]
    noup_datetimes = pd.to_datetime(sorted([inner for outer in noup_datetimes for inner in outer]))

    return up_datetimes, noup_datetimes


def season_mapping(season_str):
    if season_str == 'DJF':
        season = 'Winter'
    elif season_str == 'MAM':
        season = 'Spring'
    elif season_str == 'JJA':
        season = 'Summer'
    elif season_str == 'SON':
        season = 'Fall'
    else:
        raise ValueError(f'Invalid season string provided: {season_str}')

    return season


def setup_logger(name, loglevel, logfile):
    logger = logging.getLogger(name)

    # if the logger doesn't already exist, set it up
    if not logger.handlers:
        log_format = logging.Formatter('%(asctime)s%(module)s:%(levelname)s:%(message)s [line %(lineno)d]')
        handler = logging.FileHandler(logfile)
        handler.setFormatter(log_format)

        log_level = getattr(logging, loglevel)
        logger.setLevel(log_level)
        logger.addHandler(handler)

    return logger


def subset_wrf_grid(extent, data):
    lonx = data.XLONG
    laty = data.XLAT

    lon_idx = np.logical_and(lonx > extent[0], lonx < extent[1])
    lat_idx = np.logical_and(laty > extent[2], laty < extent[3])

    # find i and j indices of lon/lat in boundaries
    ind = np.where(np.logical_and(lat_idx, lon_idx))

    # subset data from min i,j corner to max i,j corner
    # there will be some points outside of defined boundaries because grid is not rectangular
    data_sub = np.squeeze(data)[range(np.min(ind[0]), np.max(ind[0]) + 1), range(np.min(ind[1]), np.max(ind[1]) + 1)]

    return data_sub, data_sub.XLONG, data_sub.XLAT


def wind_uv_to_dir(u, v):
    """
    Calculates the wind direction from the u and v component of wind.
    Takes into account the wind direction coordinates is different than the
    trig unit circle coordinate. If the wind direction is 360 then returns zero
    (by %360). Returns direction wind is coming FROM.
    Inputs:
    u = west/east direction (wind from the west is positive, from the east is negative)
    v = south/noth direction (wind from the south is positive, from the north is negative)
    """
    wdir = (270-np.rad2deg(np.arctan2(v, u))) % 360
    return wdir


def wind_dir_to_quadrant(wdir):
    """
    Determine directional quadrant the wind is coming from.
    Takes wind direction in degrees (0-360) and outputs quadrant:
        NE: 0-90 degrees, SE: 90-180 degrees, SW: 180-270 degrees, NW: 270-360 degrees
    Includes minimum limit (degree) of quadrant up to but not including maximum limit
    """
    wquadrant = np.array(['NA']*len(wdir))
    wquadrant[np.logical_and(wdir>=0,wdir<90)]='NE'
    wquadrant[np.logical_and(wdir>=90,wdir<180)]='SE'
    wquadrant[np.logical_and(wdir>=180,wdir<270)]='SW'
    wquadrant[np.logical_and(wdir>=270,wdir<360)]='NW'

    return wquadrant


def get_predominant_quadrant(wquadrant):
    """
    Determine predominant quadrant the wind is coming from.
    Takes array of wind quadrants and returns the one most frequently seen, plus percent of time it was observed.
    """
    quadrant_counts = wquadrant.value_counts()
    max_index = quadrant_counts.argmax()
    predominant_quadrant = quadrant_counts.keys()[max_index]
    predominant_quadrant_percent = quadrant_counts[max_index]/len(wquadrant)*100

    return predominant_quadrant, predominant_quadrant_percent


def wind_uv_to_spd(u, v):
    """
    Calculates the wind speed from the u and v wind components
    :param u: west/east direction (wind from the west is positive, from the east is negative)
    :param v: south/noth direction (wind from the south is positive, from the north is negative)
    :returns WSPD: wind speed calculated from the u and v wind components
    """
    wspd = np.sqrt(np.square(u) + np.square(v))

    return wspd


def assign_upwelling(data, upwelling_file, satellite, offset_day=True):
    """
    Assigns satellite-based upwelling to pandas dataframe
    data: pandas dataframe with 'time' column
    upwelling_file: name of csv file containing upwelling 
        (columns: Month, Day, Year, and column matching 'satellite' argument with 1=upwelling 0=no upwelling)
    satellite: name of column in upwelling_file that defines upwelling yes/no
    offset_day: whether to add one day to the dates in upwelling_file before matching to data file (default: True)
        accounts for fact that upwelling is determined by looking at maps from the day before that SST is used in the model
    """

    clms = list(data.columns)
    upwelling = pd.read_csv(upwelling_file)
    upwelling['date'] = pd.to_datetime(upwelling['Year'].astype(str)+'-'+upwelling['Month']+'-'+upwelling['Day'].astype(str))
    if offset_day:
        upwelling['date'] = upwelling['date'] + pd.Timedelta(days=1)
    data['date'] = pd.to_datetime(data['time'].dt.date)
    upwelling = upwelling[['date',satellite]]
    data = pd.merge(data, upwelling, on='date', how='left')
    data['upwelling'] = 'upwelling NA'
    data['upwelling'][data[satellite]==1] = 'upwelling present'
    data['upwelling'][data[satellite]==0] = 'upwelling absent'
    clms.append('upwelling')
    data=data[clms]

    return data


def assign_seabreeze_days(data, seabreeze_file):
    """
    Assigns seabreeze days to pandas dataframe (if seabreeze Y, ENTIRE day is defined as seabreeze)
    data: pandas dataframe with 'time' column
    seabreeze_file: name of csv file containing seabreeze dates 
        (columns: date (yyyy-mm-dd) and sea-breeze (Y=seabreeze, N=no seabreeze))
    """

    clms = list(data.columns)
    seabreeze = pd.read_csv(seabreeze_file)
    seabreeze['date'] = pd.to_datetime(seabreeze['date'], format='%Y-%m-%d')
    data['date'] = pd.to_datetime(data['time'].dt.date)
    seabreeze = seabreeze[['date','sea_breeze']]
    data = pd.merge(data, seabreeze, on='date', how='left')
    data['seabreeze'] = 'seabreeze NA'
    data['seabreeze'][data['sea_breeze']=='Y'] = 'seabreeze present'
    data['seabreeze'][data['sea_breeze']=='N'] = 'seabreeze absent'
    clms.append('seabreeze')
    data=data[clms]

    return data
