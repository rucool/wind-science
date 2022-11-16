#!/usr/bin/env python

"""
Author: Laura Nazzaro on 11/10/2022 (modified from lgarzio timeseries_plots.py)
Last modified: Lori Garzio 11/11/2022
Generate grouped summary statistics and boxplots of windspeed, power, and capacity factor
"""

import os
import glob
import numpy as np
import pandas as pd
import pickle
import matplotlib as mpl
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


def main(pfiledir, pfilepattern, out_file_dir):
    files = glob.glob(os.path.join(pfiledir, pfilepattern))
    heat_cmap = mpl.colormaps['Reds']
    heat_cmap.set_under('white')
    for f in files:
        print(f'reading {f}')
        with open(f, 'rb') as handle:
            data = pickle.load(handle)

        lease_code = data['lease'].split(' - ')[0]

        # calculate wind speed
        data['speed'] = wind_uv_to_spd(data['u'], data['v'])
        lon = np.round(np.unique(data['wrf_lon'])[0], 2)
        lat = np.round(np.unique(data['wrf_lat'])[0], 2)

        # calculate wind power and capacity factor
        power_curve = os.path.join(pfiledir,'wrf_lw15mw_power.csv')
        pc = pd.read_csv(power_curve)
        data['power'] = np.interp(data['speed'], pc['Wind Speed'], pc['Power'])
        data['capacity_factor'] = data['power'] / 15000

        # define units
        unit_labels = dict(speed='m/s',power='kW',capacity_factor=' ')
        season_order = dict(winter=1,spring=2,summer=3,fall=4)

        # set up save directories
        csv_savedir = os.path.join(out_file_dir, lease_code, 'csv')
        boxplot_savedir = os.path.join(out_file_dir, lease_code, 'boxplot')
        heatmap_savedir = os.path.join(out_file_dir, lease_code, 'heatmap')
        os.makedirs(csv_savedir, exist_ok=True)
        os.makedirs(boxplot_savedir, exist_ok=True)
        os.makedirs(heatmap_savedir, exist_ok=True)

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
        df.to_csv(os.path.join(csv_savedir, f'{lease_code}-ruwrf-timeseries_{t0.strftime("%Y%m%d")}-{t1.strftime("%Y%m%d")}.csv'))
        
        # reformat to add day of year, change to december-november year
        df['doy'] = df['time'].dt.dayofyear
        df.loc[df['month']==12, 'doy'] -= 365
        df.loc[df['month']==12, 'year'] += 1
        doy_ticks = pd.DataFrame(np.append(12,range(1,12)),columns=['month'])
        doy_ticks['day'] = 15
        doy_ticks['year'] = 2011
        doy_ticks['time'] = pd.to_datetime(doy_ticks)
        doy_ticks['doy'] = doy_ticks['time'].dt.dayofyear
        doy_ticks.loc[doy_ticks['month']==12, 'doy'] -= 365
        doy_ticks['labels']=['Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']
        t0=min(df['time'])
        t1=max(df['time'])

        keyvars = ['speed','power','capacity_factor']
        groupvars = ['year','month','season']

        for gv in groupvars:
            print(f'getting summary stats for {lease_code} by {gv}')
            all_vars = keyvars.copy()
            all_vars.append(gv)
            stats = df[all_vars].groupby(gv).describe()
            stats.to_csv(os.path.join(csv_savedir, f'{lease_code}-summary_stats_by_{gv}.csv'))
            dfgroup=df.copy()
            dfgroup['order_var'] = np.nan
            if gv in ['month','year']:
                dfgroup['order_var'] = dfgroup[gv]
                if gv=='month':
                    dfgroup.loc[dfgroup['month']==12,'order_var']=0
                    dfgroup['order_var']+=1
                    xt = list(range(1,13))
                    xl = ['Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']
                else:
                    xl = list(range(np.min(dfgroup[gv]),np.max(dfgroup[gv]+1)))
                    xt = [x+1-min(xl) for x in xl]
            elif gv=='season':
                for s in list(season_order.keys()):
                    dfgroup.loc[dfgroup['season']==s,'order_var'] = season_order[s]
                xt = list(range(1,5))
                xl = list(season_order.keys())
            bp = dfgroup.boxplot(column=keyvars, by='order_var', sharex=True, sharey=False, notch=True, showmeans=True,
                    fontsize=8, figsize=(12,4), layout=(1,3), 
                    boxprops=dict(linewidth=1.25),
                    medianprops=dict(linewidth=5,color='black'),
                    flierprops=dict(marker='x'),
                    whiskerprops=dict(linewidth=1.25),
                    meanprops=dict(marker='^', markeredgecolor='black', markerfacecolor='black'))
            for ax in bp:
                ax.set_xticks(xt,xl)
                ax.set_xlabel(None)
            plt.suptitle(None)
            plt.savefig(os.path.join(boxplot_savedir, f'{lease_code}-wind_by_{gv}_{t0.strftime("%Y%m%d")}-{t1.strftime("%Y%m%d")}.png'), dpi=300)
            plt.close()
            for kv in keyvars:
                print(f'plotting {kv} for {lease_code} by {gv}')
                fig, ax = plt.subplots(figsize=(12,7))
                dfgroup.boxplot(column=kv, by='order_var', ax=ax, notch=True, showmeans=True,
                        boxprops=dict(linewidth=1.25),
                        whiskerprops=dict(linewidth=1.25),
                        medianprops=dict(linewidth=5,color='black'),
                        flierprops=dict(marker='x'),
                        meanprops=dict(marker='^', markeredgecolor='black', markerfacecolor='black', markersize=15))
                plt.xticks(xt,xl)
                plt.xlabel(None)
                plt.ylabel(unit_labels[kv])
                plt.suptitle(f'{data["lease"]}')
                plt.title(f'{kv} {t0.strftime("%Y-%m-%d")} to {t1.strftime("%Y-%m-%d")}')
                plt.savefig(os.path.join(boxplot_savedir, f'{lease_code}-{kv}_by_{gv}_{t0.strftime("%Y%m%d")}-{t1.strftime("%Y%m%d")}.png'), dpi=300)
                plt.close()
        
        for kv in keyvars:
            if kv=='speed':
                vm = 15
            else:
                vm=8
            print(f'plotting {kv} heatmap for {lease_code}')
            fig, ax = plt.subplots(figsize=(12,7))
            plt.hexbin(dfgroup['doy'], dfgroup[kv], cmap=heat_cmap, 
                    gridsize=52, vmin=1, vmax=vm)
            plt.xticks(doy_ticks['doy'],doy_ticks['labels'])
            plt.xlabel(None)
            plt.ylabel(unit_labels[kv])
            plt.suptitle(f'{data["lease"]}')
            plt.title(f'{kv} {t0.strftime("%Y-%m-%d")} to {t1.strftime("%Y-%m-%d")}')
            plt.savefig(os.path.join(heatmap_savedir, f'{lease_code}-{kv}_heatmap_{t0.strftime("%Y%m%d")}-{t1.strftime("%Y%m%d")}.png'), dpi=300)
            plt.close()


if __name__ == '__main__':
    pckl_file_dir = '/Users/nazzaro/Documents/GitHub/wind-science/capacity_factor_analysis/files'
    pckl_file_pattern = '*.pickle'
    out_file_dir = '/Users/nazzaro/Desktop/wind_capacity'
    main(pckl_file_dir, pckl_file_pattern, out_file_dir)
