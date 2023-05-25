#! /usr/bin/env python

import logging
import numpy as np
import pandas as pd


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


def assign_upwelling(data, upwelling_file, satellite):
    """
    Assigns satellite-based upwelling to pandas dataframe
    data: pandas dataframe with 'time' column
    upwelling_file: name of csv file containing upwelling 
        (columns: Month, Day, Year, and column matching 'satellite' argument with 1=upwelling 0=no upwelling)
    satellite: name of column in upwelling_file that defines upwelling yes/no
    """

    upwelling = pd.read_csv(upwelling_file)
    upwelling['date'] = pd.to_datetime(upwelling['Year'].astype(str)+'-'+upwelling['Month']+'-'+upwelling['Day'].astype(str))
    data['date'] = pd.to_datetime(data['time'].dt.date)
    upwelling = upwelling[['date',satellite]]
    data = pd.merge(data, upwelling, on='date', how='left')
    data['upwelling'] = 'upwelling NA'
    data['upwelling'][data[satellite]==1] = 'upwelling present'
    data['upwelling'][data[satellite]==0] = 'upwelling absent'
    data=data.drop(columns=['date',satellite])

    return data