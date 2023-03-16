#! /usr/bin/env python

import logging
import numpy as np
import pandas as pd


def daterange_interval(interval, time_array):
    """
    Break up a time xarray data array into the plotting interval specified
    """
    dr_interval = []
    years = np.unique(time_array['time.year'])
    for year in years:
        year_idx = time_array['time.year'] == year
        if interval == 'monthly':
            months = np.unique(time_array['time.month'][year_idx])
            for month in months:
                month_idx = time_array[year_idx]['time.month'] == month
                time_int = time_array[year_idx][month_idx]
                dr_interval.append([np.nanmin(time_int), np.nanmax(time_int)])

        elif interval == 'seasonal':
            seasons = np.unique(time_array['time.season'][year_idx])
            for season in seasons:
                season_idx = time_array[year_idx]['time.season'] == season
                time_int = time_array[year_idx][season_idx]
                min_date = np.nanmin(time_int)
                max_date = np.nanmax(time_int)
                if season == 'DJF':  # the only season that spans 2 years
                    if np.logical_and(pd.to_datetime(min_date).month == 12, pd.to_datetime(max_date).month == 12):
                        # complete the season with the next year
                        nextyear = year + 1
                        nextyear_idx = time_array['time.year'] == nextyear
                        season_idx2 = np.logical_and(time_array[nextyear_idx]['time.season'] == season,
                                                     time_array[nextyear_idx]['time.month'] < 12)
                        time_int2 = time_array[nextyear_idx][season_idx2]
                        append_min = min_date
                        if len(time_int2) > 1:
                            max_date2 = np.nanmax(time_int2)
                            append_max = max_date2
                        else:
                            append_max = max_date
                    elif pd.to_datetime(max_date).month < 12:
                        # complete the season with the previous year
                        prevyear = year - 1
                        prevyear_idx = time_array['time.year'] == prevyear
                        season_idx2 = np.logical_and(time_array[prevyear_idx]['time.season'] == season,
                                                     time_array[prevyear_idx]['time.month'] == 12)
                        time_int2 = time_array[prevyear_idx][season_idx2]
                        append_max = max_date
                        if len(time_int2) > 0:
                            min_date2 = np.nanmin(time_int2)
                            append_min = min_date2
                        else:
                            append_min = min_date
                    else:
                        # there are 2 different seasons in this grouping, need to split them up
                        # find the earlier season
                        season1_idx = time_int['time.month'] < 12
                        season1 = time_int[season1_idx]
                        min_date = np.nanmin(season1)
                        max_date = np.nanmax(season1)

                        # complete the season with the previous year
                        prevyear = year - 1
                        prevyear_idx = time_array['time.year'] == prevyear
                        season_idx2 = np.logical_and(time_array[prevyear_idx]['time.season'] == season,
                                                     time_array[prevyear_idx]['time.month'] == 12)
                        time_int2 = time_array[prevyear_idx][season_idx2]
                        append_max = max_date
                        if len(time_int2) > 0:
                            min_date2 = np.nanmin(time_int2)
                            append_min = min_date2
                        else:
                            append_min = min_date

                        if not [append_min, append_max] in dr_interval:
                            dr_interval.append([append_min, append_max])

                        # find the later season
                        season2_idx = time_int['time.month'] == 12
                        season2 = time_int[season2_idx]
                        min_date = np.nanmin(season2)
                        max_date = np.nanmax(season2)

                        # complete the season with the next year
                        nextyear = year + 1
                        nextyear_idx = time_array['time.year'] == nextyear
                        season_idx2 = np.logical_and(time_array[nextyear_idx]['time.season'] == season,
                                                     time_array[nextyear_idx]['time.month'] < 12)
                        time_int2 = time_array[nextyear_idx][season_idx2]
                        append_min = min_date
                        if len(time_int2) > 1:
                            max_date2 = np.nanmax(time_int2)
                            append_max = max_date2
                        else:
                            append_max = max_date

                        if not [append_min, append_max] in dr_interval:
                            dr_interval.append([append_min, append_max])
                else:
                    append_min = min_date
                    append_max = max_date

                if not [append_min, append_max] in dr_interval:
                    dr_interval.append([append_min, append_max])

        else:
            raise ValueError(f'Invalid interval specified: {interval}')

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
    (by %360)
    Inputs:
    u = west/east direction (wind from the west is positive, from the east is negative)
    v = south/noth direction (wind from the south is positive, from the north is negative)
    """
    wdir = (270-np.rad2deg(np.arctan2(v, u))) % 360
    return wdir


def wind_uv_to_spd(u, v):
    """
    Calculates the wind speed from the u and v wind components
    :param u: west/east direction (wind from the west is positive, from the east is negative)
    :param v: south/noth direction (wind from the south is positive, from the north is negative)
    :returns WSPD: wind speed calculated from the u and v wind components
    """
    wspd = np.sqrt(np.square(u) + np.square(v))

    return wspd
