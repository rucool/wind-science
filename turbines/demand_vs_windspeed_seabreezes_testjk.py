#!/usr/bin/env python

"""
Author: Lori Garzio on 3/21/2022
Last modified: 3/22/2022
Plot hourly energy demand (downloaded from https://dataminer2.pjm.com/feed/hrl_load_metered) vs windspeed and estimated
wind power at 160m at the "Endurance" point in the middle of the WEA off of Atlantic City, NJ.
Seabreeze vs non-seabreeze
"""

import datetime as dt
import os
import xarray as xr
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def main(sDir, pjm_dir, sdate, edate):
    wrf = 'http://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
    point = [39.17745, -74.18033]
    plt_vars = ['windspeed', 'power']

    # for calculating power
    # power_curve = pd.read_csv('/home/lgarzio/rucool/bpu/wrf/wrf_lw15mw_power.csv')  # on server
    power_curve = pd.read_csv('/Users/garzio/Documents/rucool/bpu/wrf/wrf_lw15mw_power.csv')

    # get seabreeze dates
    seabreeze_df = pd.read_csv(os.path.join(sDir, 'radar_seabreezes_2020.csv'))
    df_sb = seabreeze_df[seabreeze_df['Seabreeze'] == 'y']
    sb_dates = np.array(pd.to_datetime(df_sb['Date']))
    sb_datetimes = [pd.date_range(pd.to_datetime(x), pd.to_datetime(x) + dt.timedelta(hours=23), freq='H') for x in
                    sb_dates]
    sb_datetimes = pd.to_datetime(sorted([inner for outer in sb_datetimes for inner in outer]))
    sb_datetimes_df = pd.DataFrame(sb_datetimes, columns=['time'])
    sb_datetimes_df['seabreeze'] = 'y'
    sb_datetimes_df.set_index('time', inplace=True)

    df_nosb = seabreeze_df[seabreeze_df['Seabreeze'] == 'n']
    nosb_dates = np.array(pd.to_datetime(df_nosb['Date']))
    nosb_datetimes = [pd.date_range(pd.to_datetime(x), pd.to_datetime(x) + dt.timedelta(hours=23), freq='H') for x
                      in nosb_dates]
    nosb_datetimes = pd.to_datetime(sorted([inner for outer in nosb_datetimes for inner in outer]))
    nosb_datetimes_df = pd.DataFrame(nosb_datetimes, columns=['time'])
    nosb_datetimes_df['seabreeze'] = 'n'
    nosb_datetimes_df.set_index('time', inplace=True)

    seabreeze_dataframe = sb_datetimes_df.append(nosb_datetimes_df)
    seabreeze_dataframe.sort_index(inplace=True)

    # get the WRF data from one point in the middle of the AC WEA at 160m
    ds = xr.open_dataset(wrf)
    ds = ds.sel(time=slice(sdate, edate))
    #ds = ds.sel(time=slice(dt.datetime(2020, 6, 1, 18, 0), dt.datetime(2020, 6, 2, 5, 0)))  # for debugging

    # calculate the sum of the absolute value distance between the model location and buoy location
    a = abs(ds['XLAT'] - point[0]) + abs(ds['XLONG'] - point[1])

    # find the indices of the minimum value in the array calculated above
    i, j = np.unravel_index(a.argmin(), a.shape)

    # get u and v component, calculate windspeed and estimated wind power
    u = ds.sel(height=160)['U'][:, i, j]
    v = ds.sel(height=160)['V'][:, i, j]
    ws = cf.wind_uv_to_spd(u, v)
    ws_df = ws.to_dataframe('windspeed')

    # estimated 15MW wind power (kW)
    power = xr.DataArray(np.interp(ws, power_curve['Wind Speed'], power_curve['Power']), coords=ws.coords)
    power_df = power.to_dataframe('power')
    power_df['power'] = power_df['power'] * .001  # convert kw to mw

    # read in the PJM energy demand data
    csv_files = glob.glob(os.path.join(pjm_dir, '*.csv'))  # '*_AE.csv'
    df = pd.DataFrame()
    for csv_file in csv_files:
        df1 = pd.read_csv(csv_file)
        df = df.append(df1)  # append to the main dataframe

    # get the total hourly demand
    df['time'] = pd.to_datetime(df['datetime_beginning_utc'])
    df.set_index('time', inplace=True)
    df_total = df.groupby('time').sum()

    merged1 = df_total.merge(power_df, on='time', how='outer')
    merged2 = merged1.merge(ws_df, on='time', how='outer')
    merged = merged2.merge(seabreeze_dataframe, on='time', how='outer')
    merged.dropna(inplace=True)

    for pv in plt_vars:
        savedir = os.path.join(sDir, f'demand_vs_{pv}')
        os.makedirs(savedir, exist_ok=True)

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 6), sharex=True, sharey=True)
        plt.subplots_adjust(left=0.08, right=0.92)

        sb = merged[merged['seabreeze'] == 'y']
        nosb = merged[merged['seabreeze'] == 'n']

        ax1.scatter(merged[pv], merged.mw, color='lightgray', s=25)
        ax2.scatter(merged[pv], merged.mw, color='lightgray', s=25)

        ax1.scatter(sb[pv], sb.mw, color='mediumblue', edgecolor='k', s=25, linewidth=.5)
        ax2.scatter(nosb[pv], nosb.mw, color='darkorange', edgecolor='k', s=25, linewidth=.5)

        ax1.set_title('Seabreeze Days: June-July 2020')
        ax2.set_title('Non-Seabreeze Days: June-July 2020')

        ax1.set_ylabel('Electricity Demand (MW)')
        if pv == 'power':
            xaxis = 'Estimated Power from 15 MW Wind Turbine (MW)'
        elif pv == 'windspeed':
            xaxis = 'Wind Speed at 160m (m/s)'
            xposition = [3, 10.9]
            for xc in xposition:
                ax1.axvline(x=xc, color='dimgray', linestyle='--')
                ax2.axvline(x=xc, color='dimgray', linestyle='--')
        ax1.set_xlabel(xaxis)
        ax2.set_xlabel(xaxis)

        save_filename = 'demand_vs_{}_seabreeze_{}-{}.png'.format(pv, sdate.strftime('%Y%m%d'), edate.strftime('%Y%m%d'))  # -AE
        sfile = os.path.join(savedir, save_filename)
        plt.savefig(sfile, dpi=300)
        plt.close()


if __name__ == '__main__':
    save_directory = '/Users/garzio/Documents/rucool/bpu/wrf/windspeed_averages'
    pjm_loc = '/Users/garzio/Documents/rucool/bpu/wrf/PJM_energy_load'
    #save_directory = '/www/home/lgarzio/public_html/bpu/windspeed_averages'  # on server
    start_date = dt.datetime(2020, 6, 1, 0, 0)  # dt.datetime(2019, 9, 1, 0, 0)
    end_date = dt.datetime(2020, 7, 31, 23, 0)  # dt.datetime(2020, 9, 1, 0, 0)
    main(save_directory, pjm_loc, start_date, end_date)
