#!/usr/bin/env python

"""
Author: Lori Garzio on 6/2/2023
Last modified: 6/2/2023
Creates power rose plots from WRF data wind farm and control runs
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


def main(args):
    start_str = args.start_str
    end_str = args.end_str
    interval = args.interval
    domain = args.domain
    save_dir = args.save_dir

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    if domain == '1km_wf2km':
        mlink = '/home/coolgroup/bpu/wrf/data/wrf_nc/1km_wf2km/1km_wf2km_160.nc'
        title_label = '1km Wind Farm'
    elif domain == '1km_ctrl':
        mlink = '/home/coolgroup/bpu/wrf/data/wrf_nc/1km_ctrl/1km_ctrl_160.nc'
        title_label = '1km Control'
    else:
        raise ValueError('Invalid domain specified')

    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=slice(start_date, end_date))

    # break up date range into the plotting interval specified
    if interval == 'none':
        intervals = [[start_date, end_date]]
    else:
        test = ds.time
        intervals = cf.daterange_interval(interval, test)
        save_dir = os.path.join(save_dir, interval)

    os.makedirs(save_dir, exist_ok=True)

    for intvl in intervals:
        sd = pd.to_datetime(intvl[0])
        ed = pd.to_datetime(intvl[1])
        dst = ds.sel(time=slice(sd, ed))
        if len(dst.time) == 0:
            raise ValueError(f'No data found for: {domain}, {sd.strftime("%Y%m%d")} to {ed.strftime("%Y%m%d")}')

        ax = new_axes()
        if interval == 'monthly':
            title_dt = sd.strftime("%b %Y")
            save_dt = sd.strftime("%Y%m%d")
        else:
            title_dt = f'{sd.strftime("%Y-%m-%d")} to {ed.strftime("%Y-%m-%d")}'
            save_dt = f'{sd.strftime("%Y%m%d")}_{ed.strftime("%Y%m%d")}'

        plt_title = f'RU-WRF {title_label}\n{title_dt}'
        sname = f'powerrose_{domain}_{save_dt}.png'
        sfile = os.path.join(save_dir, sname)

        try:
            power = dst.POWER.values / 1000  # power in kW
        except AttributeError:
            power_curve = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_lw15mw_power.csv')
            ws = cf.wind_uv_to_spd(dst.U.values, dst.V.values)
            power = np.interp(ws, power_curve['Wind Speed'], power_curve['Power'])

        power = power.flatten()
        wdir = cf.wind_uv_to_dir(dst.U.values, dst.V.values)
        wdir = wdir.flatten()

        plt_powerrose(ax, power, wdir, plt_title)

        plt.savefig(sfile, dpi=300)
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
                            default='20221031',
                            type=str,
                            help='End Date in format YYYYMMDD')

    arg_parser.add_argument('-interval',
                            dest='interval',
                            default='monthly',
                            type=str,
                            choices=['monthly', 'seasonal', 'none'],
                            help='Interval into which the time range provided is divided. If "none", the entire time '
                                 'range provided is grouped into one windrose.')

    arg_parser.add_argument('-d', '--domain',
                            dest='domain',
                            default='1km_wf2km',
                            type=str,
                            choices=['1km_wf2km', '1km_ctrl'],
                            help='Research 1km with simulated windfarm: 1km_wf2km,'
                                 'Research 1km control: 1km_ctrl')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/windturbs/plots/powerrose',
                            type=str,
                            help='Full directory path to save output plots.')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
