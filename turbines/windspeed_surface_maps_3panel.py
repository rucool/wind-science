#!/usr/bin/env python

"""
Author: Lori Garzio on 3/3/2023
Last modified: 3/3/2023
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
    heights = args.heights
    ws_lims = args.ws_lims
    save_dir = args.save_dir

    plot_turbs = '/Users/garzio/Documents/repo/rucool/wind-science/files/turbine_locations_final.csv'
    file_dir = '/home/coolgroup/ru-wrf/real-time/v4.1_parallel/processed_windturbs'

    yr = pd.to_datetime(ymd).year
    ym = ymd[0:6]

    save_dir = os.path.join(save_dir, 'surface_maps_3panel', str(yr), ym, ymd)
    os.makedirs(save_dir, exist_ok=True)

    plt_region = pf.plot_regions()
    extent = plt_region['windturb']['extent']

    files = sorted(glob.glob(os.path.join(file_dir, '1km_wf2km', ymd, '*.nc')))
    for fname in files:
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

            diff = speed - speed_ctrl
            masked_diff = np.ma.masked_inside(diff, -0.5, 0.5)

            lon = speed.XLONG.values
            lat = speed.XLAT.values

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 8), sharey=True,
                                                subplot_kw=dict(projection=ccrs.Mercator()))
            fig.suptitle(main_title, fontsize=17, y=.85)

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
            levels_ws = list(np.arange(ws_lims[0], ws_lims[1], ws_lims[2]))
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
            qs = 5
            ax1.quiver(lon[::qs, ::qs], lat[::qs, ::qs], u_std_ctrl.values[::qs, ::qs], v_std_ctrl.values[::qs, ::qs],
                       scale=30, width=.002, headlength=4, transform=ccrs.PlateCarree())
            ax2.quiver(lon[::qs, ::qs], lat[::qs, ::qs], u_std.values[::qs, ::qs], v_std.values[::qs, ::qs],
                       scale=30, width=.002, headlength=4, transform=ccrs.PlateCarree())

            # add turbine locations to wf and diff
            df = pd.read_csv(plot_turbs)
            ax2.scatter(df.lon, df.lat, s=.5, color='k', transform=ccrs.PlateCarree())
            ax3.scatter(df.lon, df.lat, s=.5, color='k', transform=ccrs.PlateCarree())

            ax1.set_title('Control')
            ax2.set_title('Wind Farm')
            ax3.set_title('Difference')
            plt.savefig(save_file, dpi=200)
            plt.close()


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description='Plot daily WRF SST input',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('ymd',
                            type=str,
                            help='Year-month-day to plot in the format YYYYmmdd (e.g. 20220101.')

    arg_parser.add_argument('-heights',
                            type=list,
                            default=[160],
                            choices=[[160], [10, 160]],
                            help='list of heights to plot')

    arg_parser.add_argument('-ws_lims',
                            type=list,
                            default=[3, 11.5, .5],
                            help='list of [minimum, maximum, interval] colorbar limits for windspeed plots')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/windturbs/plots',
                            type=str,
                            help='Full directory path to save output plots.')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
