#!/usr/bin/env python

"""
Author: Lori Garzio on 3/13/2023
Last modified: 6/20/2023
Creates wind rose plots from WRF data at user-defined time intervals, heights, and locations.
"""

import argparse
import sys
import numpy as np
import os
import xarray as xr
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from windrose import WindroseAxes
import functions.common as cf


def new_axes():
    """
    Create new wind rose axes
    """
    fig = plt.figure(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='w')
    rect = [0.15, 0.15, 0.75, 0.75]
    ax = WindroseAxes(fig, rect, facecolor='w')
    fig.add_axes(ax)
    return ax


def plt_windrose(axis, wspd, wdir, ttl):
    # set the bins for wind speeds
    b = [0, 5, 10, 15, 20, 25, 30]
    axis.bar(wdir, wspd, normed=True, bins=b, opening=1, edgecolor='black', cmap=cm.jet, nsector=36)

    # add % to y-axis labels
    newticks = ['{:.0%}'.format(x / 100) for x in axis.get_yticks()]
    axis.set_yticklabels(newticks)

    # format legend
    # move legend
    al = axis.legend(borderaxespad=-7, title='Wind Speed (m/s)')

    # replace the text in the legend
    text_str = ['0$\leq$ ws <5', '5$\leq$ ws <10', '10$\leq$ ws <15', '15$\leq$ ws <20', '20$\leq$ ws <25',
                '25$\leq$ ws <30', 'ws $\geq$30']
    for i, txt in enumerate(al.get_texts()):
        txt.set_text(text_str[i])
    plt.setp(al.get_texts(), fontsize=10)

    # add title
    axis.set_title(ttl, fontsize=14)


def main(args):
    start_str = args.start_str
    end_str = args.end_str
    interval = args.interval
    location = args.location
    height = args.height
    domain = args.domain
    save_dir = args.save_dir

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

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

    # locations of lease area centroids
    location_csv = '/home/wrfadmin/toolboxes/wind-science/files/lease_centroids.csv'
    loc_df = pd.read_csv(location_csv)
    loc_df['lease_code'] = loc_df['lease'].map(lambda x: x.split(' - ')[0].replace(' ', ''))
    df = loc_df[loc_df['lease_code'] == location]
    if len(df) == 0:
        raise ValueError('Please provide a valid lease code, found in ./files/lease_centroids.csv (i.e. "OCS-A0512")')

    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=slice(start_date, end_date))

    # break up date range into the plotting interval specified
    if interval == 'none':
        intervals = [[start_date, end_date]]
        save_dir = os.path.join(save_dir, location)
    elif interval == 'seabreeze':  # break up date range into seabreeze and non-seabreeze days
        #sb_times, nosb_times = cf.return_seabreeze_datetimes(csvfile='/Users/garzio/Documents/repo/rucool/wind-science/files/BPU_Seabreeze.csv')
        sb_times, nosb_times = cf.return_seabreeze_datetimes()
        sb_times = sb_times[np.logical_and(sb_times >= start_date, sb_times <= end_date)]
        nosb_times = nosb_times[np.logical_and(nosb_times >= start_date, nosb_times <= end_date)]
        intervals = [sb_times, nosb_times]
        save_dir = os.path.join(save_dir, location, interval)
    else:
        test = ds.time
        intervals = cf.daterange_interval(interval, test)
        save_dir = os.path.join(save_dir, location, interval)

    os.makedirs(save_dir, exist_ok=True)

    for i_intvl, intvl in enumerate(intervals):
        if interval == 'seabreeze':
            dst = ds.sel(time=intvl)
            if i_intvl == 0:
                version = 'seabreeze'
                ttl_version = 'Seabreeze Days'
            else:
                version = 'noseabreeze'
                ttl_version = 'Non-Seabreeze Days'

            if len(dst.time) == 0:
                raise ValueError(f'No data found for: {domain} {version}, {start_date.strftime("%Y%m%d")} to {end_date.strftime("%Y%m%d")}')

            # define title and save names
            title_dt = f'{ttl_version} {start_date.strftime("%Y-%m-%d")} to {end_date.strftime("%Y-%m-%d")}'
            save_dt = f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}_{version}'
        else:
            sd = pd.to_datetime(intvl[0])
            ed = pd.to_datetime(intvl[1])
            #sd = dt.datetime(2022, 3, 1, 0, 0)  # for debugging
            #ed = dt.datetime(2022, 3, 2, 0, 0)  # for debugging
            dst = ds.sel(time=slice(sd, ed))
            if len(dst.time) == 0:
                raise ValueError(f'No data found for: {domain}, {sd.strftime("%Y%m%d")} to {ed.strftime("%Y%m%d")}')

            # define title and save names
            if interval == 'monthly':
                title_dt = sd.strftime("%b %Y")
                save_dt = sd.strftime("%Y%m%d")
            else:
                title_dt = f'{sd.strftime("%Y-%m-%d")} to {ed.strftime("%Y-%m-%d")}'
                save_dt = f'{sd.strftime("%Y%m%d")}_{ed.strftime("%Y%m%d")}'

        lat = dst['XLAT']
        lon = dst['XLONG']

        if height == 10:
            u = dst['U10']
            v = dst['V10']
        else:
            u = dst.sel(height=height)['U']
            v = dst.sel(height=height)['V']

        # Find the closest model point to location
        # calculate the sum of the absolute value distance between the model location and buoy location
        a = abs(lat - df.lat.values[0]) + abs(lon - df.long.values[0])

        # find the indices of the minimum value in the array calculated above
        i, j = np.unravel_index(a.argmin(), a.shape)

        # grab the data just at that location
        usub = u[:, i, j]
        vsub = v[:, i, j]

        ax = new_axes()

        plt_title = f'RU-WRF {title_label} at {location}\n{title_dt}: {height}m'
        sname = f'windrose_{domain}_{location}_{height}m_{save_dt}.png'
        sfile = os.path.join(save_dir, sname)

        ws = cf.wind_uv_to_spd(usub.values, vsub.values)
        wdir = cf.wind_uv_to_dir(usub.values, vsub.values)

        plt_windrose(ax, ws, wdir, plt_title)

        plt.savefig(sfile, dpi=150)
        plt.close()


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-s', '--start',
                            dest='start_str',
                            default='20220301',
                            type=str,
                            help='Start Date in format YYYYMMDD')

    arg_parser.add_argument('-e', '--end',
                            dest='end_str',
                            default='20230228',
                            type=str,
                            help='End Date in format YYYYMMDD')

    arg_parser.add_argument('-interval',
                            dest='interval',
                            default='monthly',
                            type=str,
                            choices=['monthly', 'seasonal', 'none', 'seabreeze'],
                            help='Interval into which the time range provided is divided. If "none", the entire time '
                                 'range provided is grouped into one windrose. If "seabreeze" the entire time range'
                                 'provided is broken into seabreeze and non-seabreeze days.')

    arg_parser.add_argument('-location',
                            dest='location',
                            default='OCS-A0499',
                            type=str,
                            help='Lease area point at which to grab wind speeds to generate wind rose. Valid OCS lease '
                                 'codes can be found in ./files/lease_centroids.csv. Example: "OCS-A0499".')

    arg_parser.add_argument('-z', '--height',
                            dest='height',
                            default=160,
                            type=int,
                            help='Height in meters')

    arg_parser.add_argument('-d', '--domain',
                            dest='domain',
                            default='3km',
                            type=str,
                            choices=['3km', '1km_wf2km', '1km_ctrl'],
                            help='Operational: 3km, research 1km with simulated windfarm: 1km_wf2km,'
                                 'research 1km control: 1km_ctrl')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/images/windrose',
                            type=str,
                            help='Full directory path to save output plots.')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
