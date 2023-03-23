#!/usr/bin/env python

"""
Author: Laura Nazzaro on 1/12/2023 
Last modified: Laura Nazzaro 3/22/2023
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

    if len(group)>1:
        ab=False
    else:
        ab=True

    winddata = pd.DataFrame()
    powerdata = pd.DataFrame()
    for f in files:
        ds = xr.open_dataset(f)
        t0i = t0
        t1i = t1
        if not t0i:
            t0i = str(min(ds.time.data))
        if not t1i:
            t1i = str(max(ds.time.data))
        ds = ds.sel(time=slice(t0i, t1i))
        df = pd.DataFrame({'time':np.tile(ds['time'],[len(ds['points']),1]).transpose().ravel(), 'U':np.array(ds['U']).ravel(), 'V':np.array(ds['V']).ravel()})
        if 'wf' in f and 'ctrl' in f:
            print('unable to parse wind farm or control based on file name')
            return 0
        elif 'wf' in f:
            df['turbines']='wind farm'
        elif 'ctrl' in f:
            df['turbines']='control'
        else:
            print('unable to parse wind farm or control based on file name')
            return 0
        df['nturbs']=len(ds['points'])
        if winddata.empty:
            winddata = df
        else:
            winddata = pd.concat([winddata, df])
        df = pd.DataFrame({'time':ds['time']})
        if 'POWER' in ds.variables:
            df['power']=np.sum(ds['POWER'],axis=1)/1000
        else:
            u=ds['U']
            v=ds['V']
            s=wind_uv_to_spd(u,v)
            p = np.interp(s, pc['Wind Speed'], pc['Power'])
            df['power']=np.sum(p,axis=1)
        if 'wf' in f and 'ctrl' in f:
            print('unable to parse wind farm or control based on file name')
            return 0
        elif 'wf' in f:
            df['turbines']='wind farm'
        elif 'ctrl' in f:
            df['turbines']='control'
        else:
            print('unable to parse wind farm or control based on file name')
            return 0
        df['nturbs']=len(ds['points'])
        if powerdata.empty:
            powerdata = df
        else:
            powerdata = pd.concat([powerdata, df])
        ds.close()
        
    
    winddata['speed'] = wind_uv_to_spd(winddata['U'], winddata['V'])
    # data['time']=pd.to_datetime(data.index)
    winddata['year'] = winddata['time'].dt.year
    winddata['month'] = winddata['time'].dt.month
    winddata['season'] = np.nan
    winddata.loc[winddata['month']==12, 'season'] = 'winter'
    winddata.loc[winddata['month']<=2, 'season'] = 'winter'
    winddata.loc[np.logical_and(winddata['month']>=3,winddata['month']<=5), 'season'] = 'spring'
    winddata.loc[np.logical_and(winddata['month']>=6,winddata['month']<=8), 'season'] = 'summer'
    winddata.loc[np.logical_and(winddata['month']>=9,winddata['month']<=11), 'season'] = 'fall'

    nturbs=np.unique(powerdata['nturbs'])
    if len(nturbs)==1:
        max_power = max(pc['Power']*nturbs)
        powerdata['capacity_factor'] = powerdata['power'] / max_power
        powerdata['year'] = powerdata['time'].dt.year
        powerdata['month'] = powerdata['time'].dt.month
        powerdata['season'] = np.nan
        powerdata.loc[powerdata['month']==12, 'season'] = 'winter'
        powerdata.loc[powerdata['month']<=2, 'season'] = 'winter'
        powerdata.loc[np.logical_and(powerdata['month']>=3,powerdata['month']<=5), 'season'] = 'spring'
        powerdata.loc[np.logical_and(powerdata['month']>=6,powerdata['month']<=8), 'season'] = 'summer'
        powerdata.loc[np.logical_and(powerdata['month']>=9,powerdata['month']<=11), 'season'] = 'fall'

    # all_leases = '_'.join(list(np.unique(data['lease_code'])))
    # all_heights = '_'.join(list(np.unique(data['height']).astype(str)))
    all_groups = '_'.join(group)

    # set up save directories
    csv_savedir = os.path.join(mainDir, 'csv')
    boxplot_savedir = os.path.join(mainDir, 'boxplot')
    heatmap_savedir = os.path.join(mainDir, 'heatmap')
    os.makedirs(csv_savedir, exist_ok=True)
    os.makedirs(boxplot_savedir, exist_ok=True)
    os.makedirs(heatmap_savedir, exist_ok=True)

    t0=min(winddata['time']).strftime('%Y%m%dT%H%M')
    t1=max(winddata['time']).strftime('%Y%m%dT%H%M')
    t0ttl=min(winddata['time']).strftime('%Y-%m-%d %H:%M')
    t1ttl=max(winddata['time']).strftime('%Y-%m-%d %H:%M')
    # lease = list(np.unique(data['lease']))[0]

    # if len(list(np.unique(data['lease_code'])))==1:
    #     ttl = f'{lease}: \n {t0ttl} to {t1ttl}'
    # else:
    ttl=f'{t0ttl} to {t1ttl}'

    if wcsv:
        df.to_csv(os.path.join(csv_savedir, f'ruwrf-wf-timeseries_{t0}-{t1}.csv'))
    
    if 'speed' in dvs:
        dv = 'speed'
        if bp:
            boxplotname = os.path.join(boxplot_savedir,f'wf-{dv}-by-{all_groups}_{t0}-{t1}')
        else:
            boxplotname = None
        if hm:
            heatmapname = os.path.join(heatmap_savedir,f'wf-{dv}_{t0}-{t1}')
        else:
            heatmapname = None
        if gf:
            csvname = os.path.join(csv_savedir, f'wf-{dv}-by-{all_groups}_{t0}-{t1}.csv')
        else:
            csvname = None
        summary.main(winddata.copy(), dv=dv, group=group, \
            boxplotFile=boxplotname, \
            heatmapFile=heatmapname, \
            csvFile=csvname, \
            ttl=ttl, \
            all_black=ab, \
            max_speed=30)
    
    if 'power' in dvs and len(nturbs)==1:
        dv = 'power'
        if bp:
            boxplotname = os.path.join(boxplot_savedir,f'wf-cumulative{dv}-by-{all_groups}_{t0}-{t1}')
        else:
            boxplotname = None
        if hm:
            heatmapname = os.path.join(heatmap_savedir,f'wf-cumulative{dv}_{t0}-{t1}')
        else:
            heatmapname = None
        if gf:
            csvname = os.path.join(csv_savedir, f'wf-cumulative{dv}-by-{all_groups}_{t0}-{t1}.csv')
        else:
            csvname = None
        summary.main(powerdata.copy(), dv=dv, group=group, max_power=max_power, \
            boxplotFile=boxplotname, \
            heatmapFile=heatmapname, \
            csvFile=csvname, \
            all_black=ab, \
            ttl=ttl)
    



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
                            default=False,
                            type=bool,
                            help='whether to plot heatmap')
    
    arg_parser.add_argument('-gf', '--grouped_csv',
                            dest='grouped_csv',
                            default=True,
                            type=bool,
                            help='whether to write grouped csv statistics')
    
    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
