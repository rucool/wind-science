#!/usr/bin/env python

"""
Author: Lori Garzio on 3/3/2023
Last modified: 6/8/2023
Creates a 3-panel plot with surface maps of instantaneous windspeed with and without turbines, and differences
"""

import argparse
import sys
import numpy as np
import pandas as pd
import os
import glob
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cool_maps.plot as cplt
import functions.common as cf
import functions.plotting as pf
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def main(args):
    ymd = args.ymd
    expt = args.expt
    heights = args.heights
    ws_min = args.ws_min
    ws_max = args.ws_max
    ws_interval = args.ws_interval
    save_dir = args.save_dir

    turb_csv_dir = '/home/wrfadmin/toolboxes/wind-science/files'
    power_curve = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_lw15mw_power.csv')
    file_dir = '/home/coolgroup/ru-wrf/real-time/v4.1_parallel/processed_windturbs'

    yr = pd.to_datetime(ymd).year
    ym = ymd[0:6]

    plt_region = pf.plot_regions()
    if expt == '1km_wf2km':
        extent = plt_region['windturb']['extent']
        turb_csv = pd.read_csv(os.path.join(turb_csv_dir, 'turbine_locations_final.csv'))
        savestr = 'surface_maps_3panel'
        qs = 5
    elif expt == '1km_wf2km_nyb':
        extent = plt_region['windturb_nyb']['extent']
        turb_csv = pd.read_csv(os.path.join(turb_csv_dir, 'turbine_locations_final_nyb.csv'))
        savestr = 'surface_maps_3panel_nyb'
        qs = 11
    else:
        raise ValueError(f'Invalid experiment provided: {expt}')

    save_dir = os.path.join(save_dir, savestr, str(yr), ym, ymd)
    os.makedirs(save_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(file_dir, expt, ymd, '*.nc')))
    for fname in files:
        if fname.split('_')[-1] == 'H000.nc':
            continue
        f = fname.split('/')[-1]

        # find the corresponding control file
        fname_ctrl = os.path.join(file_dir, '1km_ctrl', ymd, f)

        ds = xr.open_dataset(fname)
        ds_ctrl = xr.open_dataset(fname_ctrl)
        tm = pd.to_datetime(ds.Time.values[0])

        for ht in heights:
            save_name = 'windspeed_3panel_{}m_{}_H{:03d}.png'.format(ht, tm.strftime('%Y%m%d'), tm.hour)
            save_file = os.path.join(save_dir, save_name)
            main_title = f'{tm.strftime("%Y-%m-%d %H:%M")} at {ht}m'

            if ht == 10:
                u = np.squeeze(ds['U10'])
                v = np.squeeze(ds['V10'])
                uctrl = np.squeeze(ds_ctrl['U10'])
                vctrl = np.squeeze(ds_ctrl['V10'])
            else:
                u = np.squeeze(ds.sel(height=ht)['U'])
                v = np.squeeze(ds.sel(height=ht)['V'])
                uctrl = np.squeeze(ds_ctrl.sel(height=ht)['U'])
                vctrl = np.squeeze(ds_ctrl.sel(height=ht)['V'])

            # standardize the vectors so they only represent direction
            u_std = u / cf.wind_uv_to_spd(u, v)
            v_std = v / cf.wind_uv_to_spd(u, v)

            u_std_ctrl = uctrl / cf.wind_uv_to_spd(uctrl, vctrl)
            v_std_ctrl = vctrl / cf.wind_uv_to_spd(uctrl, vctrl)

            # calculate wind speed from u and v
            speed = cf.wind_uv_to_spd(u, v)
            speed_ctrl = cf.wind_uv_to_spd(uctrl, vctrl)

            lon = speed.XLONG.values
            lat = speed.XLAT.values

            # get power from windfarm file in kW
            power = ds.POWER.values / 1000

            # calculate power from control file in kW
            # find the turbine location indices
            power_ctrl = np.array([])
            for i, row in turb_csv.iterrows():
                a = abs(lat - row.lat) + abs(lon - row.lon)
                i, j = np.unravel_index(a.argmin(), a.shape)
                power_calc = np.interp(speed_ctrl[i, j].values, power_curve['Wind Speed'], power_curve['Power'])
                power_ctrl = np.append(power_ctrl, power_calc)

            # calculate total windfarm power generated in GW
            cumulative_power = np.round(np.sum(power) / 1000000, 2)
            cumulative_power_ctrl = np.round(np.sum(power_ctrl) / 1000000, 2)
            cumulative_power_diff = np.round(cumulative_power - cumulative_power_ctrl, 2)

            diff = speed - speed_ctrl
            masked_diff = np.ma.masked_inside(diff, -0.5, 0.5)

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 8), sharey=True,
                                                subplot_kw=dict(projection=ccrs.Mercator()))
            fig.suptitle(main_title, fontsize=17, y=.86)

            # args for generating maps
            kwargs = dict()
            kwargs['landcolor'] = 'none'
            kwargs['oceancolor'] = 'none'
            kwargs['ax'] = ax1
            cplt.create(extent, **kwargs)
            kwargs['tick_label_left'] = False
            kwargs['ax'] = ax2
            cplt.create(extent, **kwargs)
            kwargs['ax'] = ax3
            cplt.create(extent, **kwargs)

            # set color maps and labels
            cmap_ws = plt.get_cmap('BuPu')
            levels_ws = list(np.arange(ws_min, ws_max + ws_interval, ws_interval))
            color_label_ws = 'Wind Speed (m/s)'

            cmap_diff = plt.get_cmap('RdBu_r')
            cmap_diff.set_bad('white')
            diff_lims = [-6, 6.5, .5]  # [min, max, interval] for difference plots
            levels_diff = list(np.arange(diff_lims[0], diff_lims[1], diff_lims[2]))
            levels_diff.remove(0.0)
            cticks = list(np.arange(-6, 7, 1))
            cticks.remove(0)
            color_label_diff = 'Wind Speed Difference (m/s)'

            # plot control
            kwargs = dict()
            kwargs['cmap'] = cmap_ws
            kwargs['clab'] = color_label_ws
            kwargs['levels'] = levels_ws
            kwargs['extend'] = 'both'
            pf.plot_contourf(fig, ax1, lon, lat, speed_ctrl, **kwargs)

            # plot wind farm
            pf.plot_contourf(fig, ax2, lon, lat, speed, **kwargs)

            # plot difference
            kwargs['cmap'] = cmap_diff
            kwargs['clab'] = color_label_diff
            kwargs['levels'] = levels_diff
            kwargs['cbar_ticks'] = cticks
            pf.plot_contourf(fig, ax3, lon, lat, masked_diff, **kwargs)

            # add WEA outline
            lease = glob.glob('/home/coolgroup/bpu/mapdata/shapefiles/BOEM-Renewable-Energy-Shapefiles-current/Wind_Lease_Outlines*.shp')[0]
            kwargs = dict()
            kwargs['edgecolor'] = 'magenta'
            pf.map_add_boem_outlines(ax1, lease, **kwargs)
            pf.map_add_boem_outlines(ax2, lease, **kwargs)
            pf.map_add_boem_outlines(ax3, lease, **kwargs)

            # add vectors to wind speed maps
            ax1.quiver(lon[::qs, ::qs], lat[::qs, ::qs], u_std_ctrl.values[::qs, ::qs], v_std_ctrl.values[::qs, ::qs],
                       scale=30, width=.002, headlength=4, transform=ccrs.PlateCarree())
            ax2.quiver(lon[::qs, ::qs], lat[::qs, ::qs], u_std.values[::qs, ::qs], v_std.values[::qs, ::qs],
                       scale=30, width=.002, headlength=4, transform=ccrs.PlateCarree())

            # add turbine locations to wf and diff
            ax2.scatter(turb_csv.lon, turb_csv.lat, s=.5, color='k', transform=ccrs.PlateCarree())
            ax3.scatter(turb_csv.lon, turb_csv.lat, s=.5, color='k', transform=ccrs.PlateCarree())

            ax1.set_title(f'Control (Power: {cumulative_power_ctrl} GW)', y=1.02)
            ax2.set_title(f'Wind Farm (Power: {cumulative_power} GW)', y=1.02)
            ax3.set_title(f'Wind Speed Difference', y=1.02)
            plt.savefig(save_file, dpi=200)
            plt.close()


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description='Plot daily WRF SST input',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('ymd',
                            type=str,
                            help='Year-month-day to plot in the format YYYYmmdd (e.g. 20220101.')

    arg_parser.add_argument('-expt',
                            type=str,
                            default='1km_wf2km',
                            choices=['1km_wf2km', '1km_wf2km_nyb'],
                            help='Experiment version to plot against the control')

    arg_parser.add_argument('-heights',
                            type=list,
                            default=[160],
                            choices=[[160], [10, 160]],
                            help='list of heights to plot')

    arg_parser.add_argument('-ws_min',
                            type=float,
                            default=3,
                            help='minimum colorbar limit for windspeed plots')

    arg_parser.add_argument('-ws_max',
                            type=float,
                            default=11,
                            help='maximum colorbar limit for windspeed plots')

    arg_parser.add_argument('-ws_interval',
                            type=float,
                            default=.5,
                            help='interval for colorbar for windspeed plots')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/windturbs/plots',
                            type=str,
                            help='Full directory path to save output plots.')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
