#!/usr/bin/env python

"""
Author: Laura Nazzaro on 11/10/2022 (modified from lgarzio timeseries_plots.py)
Last modified: Laura Nazzaro 11/10/2022
Generate grouped summary statistics and boxplots of windspeed, power, and capacity factor
"""

import os
import glob
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console


def wind_uv_to_spd(u, v):
    """
    Calculates the wind speed from the u and v wind components
    :param u: west/east direction (wind from the west is positive, from the east is negative)
    :param v: south/noth direction (wind from the south is positive, from the north is negative)
    :returns wspd: wind speed calculated from the u and v wind components
    """
    wspd = np.sqrt(np.square(u) + np.square(v))

    return wspd


def main(pfiledir, out_file_dir):
    files = glob.glob(os.path.join(pfiledir, '*.pickle'))
    for f in files:
        print(f'reading {f}')
        with open(f, 'rb') as handle:
            data = pickle.load(handle)

        # calculate wind speed
        data['speed'] = wind_uv_to_spd(data['u'], data['v'])
        lon = np.round(np.unique(data['wrf_lon'])[0], 2)
        lat = np.round(np.unique(data['wrf_lat'])[0], 2)

        # calculate wind power and capacity factor
        power_curve = os.path.join(pfiledir, 'wrf_lw15mw_power.csv')
        pc = pd.read_csv(power_curve)
        data['power'] = np.interp(data['speed'], pc['Wind Speed'], pc['Power'])
        data['capacity_factor'] = data['power'] / 15000

        # assign to dataframe and add month and season
        df = pd.DataFrame.from_dict(data)
        df['year'] = df['time'].dt.year
        df['month'] = df['time'].dt.month
        df['season'] = np.nan
        df.loc[df['month']==12, 'season'] = 'winter'
        df.loc[df['month']<=2, 'season'] = 'winter'
        df.loc[np.logical_and(df['month']>=3,df['month']<=5), 'season'] = 'spring'
        df.loc[np.logical_and(df['month']>=6,df['month']<=8), 'season'] = 'summer'
        df.loc[np.logical_and(df['month']>=9,df['month']<=11), 'season'] = 'fall'
        df.to_csv(os.path.join(out_file_dir,'csv',f'{data["code"]}-ruwrf-timeseries.csv'))

        keyvars = ['speed','power','capacity_factor']
        groupvars = ['month','season']

        for gv in groupvars:
            print(f'getting summary stats for {data["code"]} by {gv}')
            all_vars = keyvars.copy()
            all_vars.append(gv)
            stats = df[all_vars].groupby(gv).describe()
            stats.to_csv(os.path.join(out_file_dir,'csv',f'{data["code"]}-summary_stats_by_{gv}.csv'))
            #fig, ax = plt.subplots(figsize=(12,7))
            df.boxplot(column=keyvars, by=gv, sharex=True, sharey=False, notch=True, fontsize=8, figsize=(12,4), layout=(1,3), boxprops=dict(linewidth=1.25),medianprops=dict(linewidth=3,color='green'),flierprops=dict(marker='x'),whiskerprops=dict(linewidth=1.25))
            plt.suptitle('') 
            plt.savefig(os.path.join(out_file_dir,'images',f'{data["code"]}-wind_by_{gv}.png'),dpi=300)
            plt.close()
            for kv in keyvars:
                print(f'plotting {kv} for {data["code"]} by {gv}')
                fig, ax = plt.subplots(figsize=(12,7))
                df.boxplot(column=kv, by=gv, ax=ax, notch=True, boxprops=dict(linewidth=1.25),whiskerprops=dict(linewidth=1.25),medianprops=dict(linewidth=3,color='green'),flierprops=dict(marker='x'))
                plt.suptitle('') 
                plt.savefig(os.path.join(out_file_dir,'images',f'{data["code"]}-{kv}_by_{gv}.png'),dpi=300)
                plt.close()




if __name__ == '__main__':
    pckl_file_dir = '/Users/nazzaro/Documents/GitHub/wind-science/capacity_factor_analysis/files'
    out_file_dir = '/Users/nazzaro/Desktop/wind_capacity'
    main(pckl_file_dir, out_file_dir)
