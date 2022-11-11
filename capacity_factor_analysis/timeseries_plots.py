#!/usr/bin/env python

"""
Author: Lori Garzio on 11/10/2022
Last modified: Lori Garzio 11/10/2022
Generate timeseries plots of windspeed, power, and capacity factor
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


def main(pfiledir):
    files = glob.glob(os.path.join(pfiledir, '*.pickle'))
    for f in files:
        with open(f, 'rb') as handle:
            data = pickle.load(handle)

        # calculate and plot wind speed
        speed = wind_uv_to_spd(data['u'], data['v'])
        lon = np.round(np.unique(data['wrf_lon'])[0], 2)
        lat = np.round(np.unique(data['wrf_lat'])[0], 2)
        save_file = os.path.join(pfiledir, f'{data["code"]}-windspeed_ts.png')

        fig, ax = plt.subplots(figsize=(12, 7))
        ax.plot(data['time'], speed)
        ax.set_ylabel('Windspeed (m/s)')
        ax.set_title(f'Windspeed {data["height"]}m at center of WEA {lon, lat}\n{data["company"]}')

        plt.savefig(save_file, dpi=300)
        plt.close()

        # calculate and plot wind power
        power_curve = '/Users/garzio/Documents/rucool/bpu/wrf/wrf_lw15mw_power.csv'
        pc = pd.read_csv(power_curve)
        power = np.interp(speed, pc['Wind Speed'], pc['Power'])

        fig, ax = plt.subplots(figsize=(12, 7))
        ax.plot(data['time'], power)
        ax.set_ylabel('Estimated 15 MW Wind Power (kW)')
        ax.set_title(f'Est Power {data["height"]}m at center of WEA {lon, lat}\n{data["company"]}')

        save_file = os.path.join(pfiledir, f'{data["code"]}-windpower_ts.png')
        plt.savefig(save_file, dpi=300)
        plt.close()

        capacity_factor = power / 15000

        fig, ax = plt.subplots(figsize=(12, 7))
        ax.plot(data['time'], capacity_factor)
        ax.set_ylabel('Estimated Capacity Factor')
        ax.set_title(f'Est Capacity Factor {data["height"]}m at center of WEA {lon, lat}\n{data["company"]}')

        save_file = os.path.join(pfiledir, f'{data["code"]}-capacityfactor_ts.png')
        plt.savefig(save_file, dpi=300)
        plt.close()


if __name__ == '__main__':
    pckl_file_dir = '/Users/garzio/Documents/rucool/bpu/wrf/capacity_factor_analysis'
    main(pckl_file_dir)
