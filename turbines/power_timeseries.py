#!/usr/bin/env python

"""
Author: Lori Garzio on 2/22/2023
Last modified: 3/22/2023
Plot timeseries of cumulative wind farm power (hourly and daily), average wind speed and direction from RU-WRF
simulated turbines vs control
"""

import argparse
import sys
import os
import datetime as dt
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import functions.common as cf
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def format_date_axis(axis):
    #datef = mdates.DateFormatter('%m-%d\n%H:%M')
    datef = mdates.DateFormatter('%b-%d\n%Y')
    axis.xaxis.set_major_formatter(datef)


def main(args):
    start_str = args.start
    end_str = args.end
    fwf = args.file_wf
    fctrl = args.file_ctrl
    save_dir = args.save_dir

    save_dir = os.path.join(save_dir, 'timeseries')
    os.makedirs(save_dir, exist_ok=True)
    power_curve = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_lw15mw_power.csv')

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    # subset time range
    ds = xr.open_dataset(fwf)
    ds_ctrl = xr.open_dataset(fctrl)

    ds = ds.sel(time=slice(start_date, end_date))
    ds_ctrl = ds_ctrl.sel(time=slice(start_date, end_date))

    # get cumulative power for the entire windfarm at each hour
    # power output by the model and power calculated using windspeeds
    cumsum_power_kw = ds.POWER.sum(dim='points') / 1000  # power in kW
    speed = cf.wind_uv_to_spd(ds.U, ds.V)
    power_calc = xr.DataArray(np.interp(speed, power_curve['Wind Speed'], power_curve['Power']), coords=speed.coords,
                              dims=speed.dims)
    cumsum_power_kw_calc = power_calc.sum(dim='points')

    # calculate average wind speed and direction
    speed_avg = speed.mean(dim='points')
    uavg = ds.U.mean(dim='points')
    vavg = ds.V.mean(dim='points')
    direction_avg = cf.wind_uv_to_dir(uavg, vavg)

    # calculate power for the windfarm location for the control run using wind speeds
    speed_ctrl = cf.wind_uv_to_spd(ds_ctrl.U, ds_ctrl.V)
    power_calc_ctrl = xr.DataArray(np.interp(speed_ctrl, power_curve['Wind Speed'], power_curve['Power']),
                                   coords=speed_ctrl.coords, dims=speed_ctrl.dims)
    cumsum_power_kw_calc_ctrl = power_calc_ctrl.sum(dim='points')

    # calculate daily cumulative power
    cumsum_power_kw_daily = cumsum_power_kw.resample(time='D', keep_attrs=True).sum()
    cumsum_power_kw_calc_ctrl_daily = cumsum_power_kw_calc_ctrl.resample(time='D', keep_attrs=True).sum()

    # calculate average wind speed and direction
    speed_avg_ctrl = speed_ctrl.mean(dim='points')
    uavg_ctrl = ds_ctrl.U.mean(dim='points')
    vavg_ctrl = ds_ctrl.V.mean(dim='points')
    direction_avg_ctrl = cf.wind_uv_to_dir(uavg_ctrl, vavg_ctrl)

    overall_sum_gw = np.round(np.sum(cumsum_power_kw.values) / 1000000, 2)
    overall_sum_calc_gw = np.round(np.sum(cumsum_power_kw_calc.values) / 1000000, 2)
    overall_sum_ctrl_gw = np.round(np.sum(cumsum_power_kw_calc_ctrl.values) / 1000000, 2)

    print(f'Total power: simulated wind farm {overall_sum_gw} GW')
    print(f'Total power: simluated wind farm calculated {overall_sum_calc_gw} GW')
    print(f'Total power: control {overall_sum_ctrl_gw} GW')

    # plot power - hourly
    fig, ax = plt.subplots(figsize=(13, 7))

    ax.plot(cumsum_power_kw_calc_ctrl.time, cumsum_power_kw_calc_ctrl / 1000000, color='#d95f02', label='Control')  # orange
    ax.plot(cumsum_power_kw.time, cumsum_power_kw / 1000000, color='#7570b3', label='Turbines')  # purple

    ax.set_ylabel('Wind Farm Power (GW)')
    ax.set_xlabel('Time (GMT)')
    ax.legend(loc='best', fontsize=10)

    format_date_axis(ax)
    ax.set_ylim([-0.1, 3.1])

    title = f'Simulated Turbines Total Power {overall_sum_gw} GW; Control Total Power {overall_sum_ctrl_gw} GW'
    plt.title(title)

    save_file = os.path.join(save_dir, f'turbine_vs_control_timeseries_power_hourly_{start_str}_{end_str}.png')
    plt.savefig(save_file, dpi=200)
    plt.close()

    # plot power - daily
    fig, ax = plt.subplots(figsize=(13, 7))

    ax.plot(cumsum_power_kw_calc_ctrl_daily.time, cumsum_power_kw_calc_ctrl_daily / 1000000, color='#d95f02', label='Control')  # orange
    ax.plot(cumsum_power_kw_daily.time, cumsum_power_kw_daily / 1000000, color='#7570b3', label='Turbines')  # purple

    ax.set_ylabel('Wind Farm Power (GW)')
    ax.set_xlabel('Time (GMT)')
    ax.legend(loc='best', fontsize=10)

    format_date_axis(ax)
    #ax.set_ylim([-0.1, 3.1])

    title = f'Simulated Turbines Total Power {overall_sum_gw} GW; Control Total Power {overall_sum_ctrl_gw} GW'
    plt.title(title)

    save_file = os.path.join(save_dir, f'turbine_vs_control_timeseries_power_daily_{start_str}_{end_str}.png')
    plt.savefig(save_file, dpi=200)
    plt.close()

    # plot power calculated vs power output by the model
    fig, ax = plt.subplots(figsize=(13, 7))

    ax.plot(cumsum_power_kw_calc.time, cumsum_power_kw_calc / 1000000, color='#d95f02', label='Turbines Calculated')  # orange
    ax.plot(cumsum_power_kw.time, cumsum_power_kw / 1000000, color='#7570b3', label='Turbines Output')  # purple

    ax.set_ylabel('Wind Farm Power (GW)')
    ax.set_xlabel('Time (GMT)')
    ax.legend(loc='best', fontsize=10)

    format_date_axis(ax)

    save_file = os.path.join(save_dir, f'turbine_timeseries_power_vs_calculated_{start_str}_{end_str}.png')
    plt.savefig(save_file, dpi=200)
    plt.close()

    # plot average windspeed for turbine locations
    fig, ax = plt.subplots(figsize=(13, 7))

    plt.axhline(y=3, ls='--', c='gray')
    plt.axhline(y=10.9, ls='--', c='gray')

    ax.plot(speed_avg_ctrl.time, speed_avg_ctrl, color='#d95f02', label='Control')  # orange
    ax.plot(speed_avg.time, speed_avg, color='#7570b3', label='Turbines')  # purple

    ax.set_ylabel('Wind Farm Average Wind Speed (m/s)')
    ax.set_xlabel('Time (GMT)')
    ax.legend(loc='best', fontsize=10)
    #ax.set_ylim(0, 20)

    format_date_axis(ax)

    save_file = os.path.join(save_dir, f'turbine_vs_control_timeseries_windspeed_{start_str}_{end_str}.png')
    plt.savefig(save_file, dpi=200)
    plt.close()

    # plot average wind direction for turbine locations
    fig, ax = plt.subplots(figsize=(13, 7))

    ax.plot(direction_avg_ctrl.time, direction_avg_ctrl, color='#d95f02', label='Control')  # orange
    ax.plot(direction_avg.time, direction_avg, color='#7570b3', label='Turbines')  # purple

    ax.set_ylabel('Wind Farm Average Wind Direction (Degrees)')
    ax.set_xlabel('Time (GMT)')
    ax.legend(loc='best', fontsize=10)

    format_date_axis(ax)

    save_file = os.path.join(save_dir, f'turbine_vs_control_timeseries_winddirection_{start_str}_{end_str}.png')
    plt.savefig(save_file, dpi=200)
    plt.close()

    # plot average wind direction as a quiver plot for turbine locations
    fig, (ax1, ax2) = plt.subplots(2, figsize=(13, 7), sharex=True, sharey=True)

    #ax.quiver(direction_avg_ctrl.time, 0, uavg_ctrl, vavg_ctrl, headlength=1, headaxislength=1, width=0.004, alpha=0.5)
    ax1.quiver(direction_avg_ctrl.time, 0, uavg_ctrl, vavg_ctrl, scale=120, alpha=0.5, label='Control')
    ax2.quiver(direction_avg.time, 0, uavg, vavg, scale=120, alpha=0.5, label='Turbines')

    ax2.set_xlabel('Time (GMT)')
    ax1.set_title('Average Wind Direction: Control')
    ax2.set_title('Average Wind Direction: Wind Farm')

    format_date_axis(ax2)

    save_file = os.path.join(save_dir, f'quiver_timeseries_{start_str}_{end_str}.png')
    plt.savefig(save_file, dpi=200)
    plt.close()


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-s', '--start',
                            dest='start',
                            default='20220801',
                            type=str,
                            help='Start Date in format YYYYMMDD ')

    arg_parser.add_argument('-e', '--end',
                            dest='end',
                            default='20220831',
                            type=str,
                            help='End Date in format YYYYMMDD')

    arg_parser.add_argument('-file_wf',
                            dest='file_wf',
                            default='/home/coolgroup/bpu/wrf/data/wrf_nc/1km_wf2km/1km_wf2km_160.nc',
                            type=str,
                            help='Path to file containing time-series subset of RU-WRF 1km run with simulated turbines.')

    arg_parser.add_argument('-file_ctrl',
                            dest='file_ctrl',
                            default='/home/coolgroup/bpu/wrf/data/wrf_nc/1km_ctrl/1km_ctrl_160.nc',
                            type=str,
                            help='Domain: 1km with simulated windfarm (1km_wf2km), 1km without simulated windfarm (1km_ctrl)')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/windturbs/plots',
                            type=str,
                            help='Full directory path to save output files.')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))

