#!/usr/bin/env python

"""
Author: Laura Nazzaro on 1/12/2023 
Last modified: Laura Nazzaro 1/12/2023
Generate grouped summary statistics and boxplots of windspeed, power, and capacity factor
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from functions import grouped_data_patterns as summary


def wind_uv_to_spd(u, v):
    """
    Calculates the wind speed from the u and v wind components
    :param u: west/east direction (wind from the west is positive, from the east is negative)
    :param v: south/noth direction (wind from the south is positive, from the north is negative)
    :returns wspd: wind speed calculated from the u and v wind components
    """
    wspd = np.sqrt(np.square(u) + np.square(v))

    return wspd

def main(args):
    t1 = args.end
    t0 = args.start
    mainDir = args.main_dir
    files = args.files.split(',')
    group = args.groups.split(',')
    wcsv = args.write_csv
    pc = pd.read_csv(args.power_curve)
    dvs = args.response_variable.split(',')
    bp = args.plot_boxplot
    hm = args.plot_heatmap
    gf = args.grouped_csv

    data = pd.DataFrame()
    for f in files:
        ds = xr.open_dataset(f)
        t0i = t0
        t1i = t1
        if not t0i:
            t0i = str(min(ds.time.data))
        if not t1i:
            t1i = str(max(ds.time.data))
        ds = ds.sel(time=slice(t0i, t1i))
        df = ds.to_dataframe()
        df['lat'] = ds.point_lat
        df['lon'] = ds.point_lon
        df['company'] = ds.company
        df['lease'] = ds.lease
        df['lease_code'] = ds.lease_code
        df['state'] = ds.state
        if data.empty:
            data = df
        else:
            data = pd.concat([data, df])
    
    max_power = max(pc['Power'])
    data['speed'] = wind_uv_to_spd(data['U'], data['V'])
    data['power'] = np.interp(data['speed'], pc['Wind Speed'], pc['Power'])
    data['capacity_factor'] = data['power'] / max_power
    data['time']=pd.to_datetime(data.index)
    data['year'] = data['time'].dt.year
    data['month'] = data['time'].dt.month
    data['season'] = np.nan
    data.loc[data['month']==12, 'season'] = 'winter'
    data.loc[data['month']<=2, 'season'] = 'winter'
    data.loc[np.logical_and(data['month']>=3,data['month']<=5), 'season'] = 'spring'
    data.loc[np.logical_and(data['month']>=6,data['month']<=8), 'season'] = 'summer'
    data.loc[np.logical_and(data['month']>=9,data['month']<=11), 'season'] = 'fall'

    all_leases = '_'.join(list(np.unique(data['lease_code'])))
    all_heights = '_'.join(list(np.unique(data['height']).astype(str)))
    all_groups = '_'.join(group)

    # set up save directories
    csv_savedir = os.path.join(mainDir, 'csv')
    boxplot_savedir = os.path.join(mainDir, 'boxplot')
    heatmap_savedir = os.path.join(mainDir, 'heatmap')
    os.makedirs(csv_savedir, exist_ok=True)
    os.makedirs(boxplot_savedir, exist_ok=True)
    os.makedirs(heatmap_savedir, exist_ok=True)

    t0=min(data['time']).strftime('%Y%m%dT%H%M')
    t1=max(data['time']).strftime('%Y%m%dT%H%M')

    if wcsv:
        df.to_csv(os.path.join(csv_savedir, f'{all_leases}-ruwrf-timeseries_{all_heights}_{t0}-{t1}.csv'))
    
    for dv in dvs:
        if bp:
            boxplotname = os.path.join(boxplot_savedir,f'{all_leases}-{dv}-by-{all_groups}_{all_heights}_{t0}-{t1}')
        else:
            boxplotname = None
        if hm:
            heatmapname = os.path.join(heatmap_savedir,f'{all_leases}-{dv}_{all_heights}_{t0}-{t1}')
        else:
            heatmapname = None
        if gf:
            csvname = os.path.join(csv_savedir, f'{all_leases}-{dv}-by-{all_groups}_{all_heights}_{t0}-{t1}.csv')
        else:
            csvname = None
        summary.main(data.copy(), dv=dv, group=group, \
            boxplotFile=boxplotname, \
            heatmapFile=heatmapname, \
            csvFile=csvname)
    



if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-e', '--end',
                            dest='end',
                            default=None,
                            type=str,
                            help='End Date in format YYYYMMDDTHH:MM')
    
    arg_parser.add_argument('-s', '--start',
                            dest='start',
                            default=None,
                            type=str,
                            help='Start Date in format YYYYMMDDTHH:MM')
    
    arg_parser.add_argument('-d', '--main_dir',
                            dest='main_dir',
                            default=None,
                            type=str,
                            help='Directory to add csv files and images to')

    arg_parser.add_argument('-f', '--files',
                            dest='files',
                            default=None,
                            type=str,
                            help='List of files to grab data from, comma-separated with no spaces')

    arg_parser.add_argument('-w', '--write_csv',
                            dest='write_csv',
                            default=False,
                            type=bool,
                            help='Whether to write csv file with all data used')
    
    arg_parser.add_argument('-p', '--power_curve',
                            dest='power_curve',
                            default=None,
                            type=str,
                            help='File containing power curve')
    
    arg_parser.add_argument('-g', '--groups',
                            dest='groups',
                            default=None,
                            type=str,
                            help='variable(s) to group by - if multiple, separated by commas (no spaces)')
    
    arg_parser.add_argument('-v', '--response_variable',
                            dest='response_variable',
                            default='power,speed',
                            type=str,
                            help='response variable(s) to use (speed, power, or capacity_factor) - if multiple, separated by commas (no spaces)')
    
    arg_parser.add_argument('-bp', '--boxplot',
                            dest='plot_boxplot',
                            default=True,
                            type=bool,
                            help='whether to plot boxplot')
    
    arg_parser.add_argument('-hm', '--heatmap',
                            dest='plot_heatmap',
                            default=True,
                            type=bool,
                            help='whether to plot heatmap')
    
    arg_parser.add_argument('-gf', '--grouped_csv',
                            dest='grouped_csv',
                            default=True,
                            type=bool,
                            help='whether to write grouped csv statistics')
    
    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
