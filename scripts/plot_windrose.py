#!/usr/bin/env python

"""
Author: Lori Garzio on 3/13/2023
Last modified: 3/26/2024
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
    # rect = [0.15, 0.15, 0.75, 0.75]
    rect = [0.15, 0.13, 0.72, 0.75]
    ax = WindroseAxes(fig, rect, facecolor='w')
    fig.add_axes(ax)
    return ax


def plt_powerrose(axis, power, wdir, ttl):
    # set the bins
    b = [0, 2500, 5000, 7500, 10000, 12500, 15000]
    axis.bar(wdir, power, normed=True, bins=b, opening=1, edgecolor='black', cmap=cm.jet, nsector=36)

    # add % to y-axis labels
    newticks = ['{:.0%}'.format(x / 100) for x in axis.get_yticks()]
    axis.set_yticklabels(newticks)

    # format legend
    # move legend
    al = axis.legend(borderaxespad=-7, title='Power (kW)')

    # replace the text in the legend
    text_str = ['0$\leq$ p <2500', '2,500$\leq$ p <5000', '5,000$\leq$ p <7500', '7,500$\leq$ p <10,000',
                '10,000$\leq$ p <12,500', '12,500$\leq$ p <15,000', 'p $\geq$15,000']
    for i, txt in enumerate(al.get_texts()):
        txt.set_text(text_str[i])
    plt.setp(al.get_texts(), fontsize=10)

    # add title
    axis.set_title(ttl, fontsize=14)


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
    plot_power = args.plot_power
    subset_dir = args.subset_dir
    save_dir = args.save_dir

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    if domain == '3km':
        if subset_dir:
            mlink = f'/home/coolgroup/bpu/wrf/data/wrf_nc/wea_centroids/3km/{location}_160_3km.nc'
            height = 160  # height must be 160m if using subset files
        else:
            mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
        title_label = '3km'
    elif domain == '1km_wf2km':
        if subset_dir:
            ValueError('subset_dir must be set to False when plotting 1km_wf2km data')
        else:
            mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_wf2km_processed/WRF_4.1_1km_with_Wind_Farm_Processed_Dataset_Best'
        title_label = '1km Wind Farm'
    elif domain == '1km_ctrl':
        if subset_dir:
            ValueError('subset_dir must be set to False when plotting 1km_ctrl data')
        else:
            mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_ctrl_processed/WRF_4.1_1km_Control_Processed_Dataset_Best'
        title_label = '1km Control'
    else:
        raise ValueError('Invalid domain specified')

    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=slice(start_date, end_date))

    # locations of lease area centroids
    location_csv = '/home/wrfadmin/toolboxes/wind-science/files/lease_centroids.csv'
    loc_df = pd.read_csv(location_csv)
    loc_df['lease_code'] = loc_df['lease'].map(lambda x: x.split(' - ')[0].replace(' ', ''))
    df = loc_df[loc_df['lease_code'] == location]
    if len(df) == 0:
        raise ValueError('Please provide a valid lease code, found in ./files/lease_centroids.csv (i.e. "OCS-A0512")')

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
    elif interval == 'upwelling':  # break up date range into upwelling and non-upwelling days
        up_times, noup_times = cf.return_upwelling_datetimes()
        up_times = up_times[np.logical_and(up_times >= start_date, up_times <= end_date)]
        noup_times = noup_times[np.logical_and(noup_times >= start_date, noup_times <= end_date)]
        intervals = [up_times, noup_times]
    elif interval == 'summers':
        months = ds['time.month']
        intervals = list(np.where((months == 6) | (months == 7) | (months == 8)))
    elif interval == 'months':
        intervals = []
        for month in np.unique(ds['time.month']):
            intervals.append(list(np.where(ds['time.month'] == month)[0]))
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
                ttl_version = f'Seabreeze Days (n={int(len(dst.time) / 24)})'
            else:
                version = 'noseabreeze'
                ttl_version = f'Non-Seabreeze Days (n={int(len(dst.time) / 24)})'

            if len(dst.time) == 0:
                raise ValueError(f'No data found for: {domain} {version}, {start_date.strftime("%Y%m%d")} to {end_date.strftime("%Y%m%d")}')

            # define title and save names
            if np.logical_and(start_str == '20200601', end_str == '20220831'):
                title_dt = f'{ttl_version}\nJune-July-Aug 2020-2022'
                save_dt = f'summers2020_2022_{version}'
            else:
                title_dt = f'{ttl_version}\n{start_date.strftime("%Y-%m-%d")} to {end_date.strftime("%Y-%m-%d")}'
                save_dt = f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}_{version}'
        elif interval == 'upwelling':
            dst = ds.sel(time=intvl)
            if i_intvl == 0:
                version = 'upwell'
                ttl_version = f'Upwelling Days (n={int(len(dst.time) / 24)})'
            else:
                version = 'noupwell'
                ttl_version = f'Non-Upwelling Days (n={int(len(dst.time) / 24)})'

            if len(dst.time) == 0:
                raise ValueError(
                    f'No data found for: {domain} {version}, {start_date.strftime("%Y%m%d")} to {end_date.strftime("%Y%m%d")}')

            # define title and save names
            if np.logical_and(start_str == '20190601', end_str == '20220831'):
                title_dt = f'{ttl_version}\nJune-July-Aug 2019-2022'
                save_dt = f'summers2019_2022_{version}'
            else:
                title_dt = f'{ttl_version}\n{start_date.strftime("%Y-%m-%d")} to {end_date.strftime("%Y-%m-%d")}'
                save_dt = f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}_{version}'
        elif interval == 'summers':
            dst = ds.isel(time=intvl)
            min_year = np.nanmin(dst['time.year'])
            max_year = np.nanmax(dst['time.year'])
            title_dt = f'June-July-Aug {str(min_year)}-{str(max_year)}'
            save_dt = f'summers{str(min_year)}_{str(max_year)}'
        elif interval == 'months':
            dst = ds.isel(time=intvl)
            years = '-'.join(str(x) for x in np.unique(dst['time.year']).tolist())
            month_num = np.unique(dst['time.month'])[0]
            month_name = pd.to_datetime(dst.time.values[0]).strftime("%b")
            title_dt = f'{month_name} {years}'
            save_dt = f'months_{str(month_num).zfill(2)}_{month_name}_{years}'
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
                title_dt = f'{sd.strftime("%b %Y")}: {int(len(dst.time)/24)} days'
                save_dt = sd.strftime("%Y%m%d")
            elif interval == 'seasonal':
                season = cf.season_mapping(np.unique(dst['time.season'])[0])
                years = '-'.join(str(x) for x in np.unique(dst['time.year']))
                title_dt = f'{season} {years}: {int(len(dst.time) / 24)} days'
                save_dt = f'{sd.strftime("%Y%m%d")}_{ed.strftime("%Y%m%d")}'
            else:
                title_dt = f'{sd.strftime("%Y-%m-%d")} to {ed.strftime("%Y-%m-%d")}: {int(len(dst.time)/24)} days'
                save_dt = f'{sd.strftime("%Y%m%d")}_{ed.strftime("%Y%m%d")}'

        if np.logical_and(subset_dir, domain == '3km'):  # if using the 3km subset file that's already specific to a centroid
            usub = dst.U
            vsub = dst.V
        else:
            lat = dst['XLAT']
            lon = dst['XLONG']

            if height == 10:
                u = dst['U10']
                v = dst['V10']
            else:
                try:
                    u = dst.sel(height=height)['U']
                    v = dst.sel(height=height)['V']
                except KeyError:  # the file only has one height (160m)
                    u = dst.U
                    v = dst.V

            # Find the closest model point to location
            # calculate the sum of the absolute value distance between the model location and buoy location
            a = abs(lat - df.lat.values[0]) + abs(lon - df.long.values[0])

            # find the indices of the minimum value in the array calculated above
            i, j = np.unravel_index(a.argmin(), a.shape)

            # grab the data at that location
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

        if plot_power:
            power_curve = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_lw15mw_power.csv')

            power = np.interp(ws, power_curve['Wind Speed'], power_curve['Power'])

            ax = new_axes()
            plt_powerrose(ax, power, wdir, plt_title)

            sname = f'powerrose_{domain}_{location}_{height}m_{save_dt}.png'
            sfile = os.path.join(save_dir, sname)
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
                            choices=['monthly', 'seasonal', 'none', 'seabreeze', 'summers', 'months', 'upwelling'],
                            help='Interval into which the time range provided is divided. If "none", the entire time '
                                 'range provided is grouped into one windrose. If "seabreeze" the entire time range'
                                 'provided is broken into seabreeze and non-seabreeze days. If "summers" the entire'
                                 'time range provided is subset for summers only (June-July-Aug). If "months" the '
                                 'entire time range provided is broken into months and multiple years are included '
                                 'in the same month plot.')

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

    arg_parser.add_argument('-plot_power',
                            default=False,
                            type=bool,
                            choices=[True, False],
                            help='Option to plot power roses, default is False')

    arg_parser.add_argument('-subset_dir',
                            default=False,
                            type=bool,
                            choices=[True, False],
                            help='Option to use a previously-generated subset file in '
                                 '/home/coolgroup/bpu/wrf/data/wrf_nc/wea_centroids/3km. If set to True, height must '
                                 'be 160m and domain must be 3km. Default is False')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/images/windrose',
                            type=str,
                            help='Full directory path to save output plots.')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
