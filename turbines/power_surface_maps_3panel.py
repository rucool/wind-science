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
    power_min = args.power_min
    power_max = args.power_max
    power_interval = args.power_interval
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
        savestr = 'power_maps_3panel'
        qs = 5
    elif expt == '1km_wf2km_nyb':
        extent = plt_region['windturb_nyb']['extent']
        turb_csv = pd.read_csv(os.path.join(turb_csv_dir, 'turbine_locations_final_nyb.csv'))
        savestr = 'power_maps_3panel_nyb'
        qs = 11
    else:
        raise ValueError(f'Invalid experiment provided: {expt}')

    save_dir = os.path.join(save_dir, savestr, str(yr), ym, ymd)
    os.makedirs(save_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(file_dir, expt, ymd, '*.nc')))
    
    # initialize totals for the 24 hour period
    cumulative_power_total = 0
    cumulative_power_ctrl_total = 0
    cumulative_power_diff_total = 0
    
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
            save_name = 'power_3panel_{}m_{}_H{:03d}.png'.format(ht, tm.strftime('%Y%m%d'), tm.hour)
            save_file = os.path.join(save_dir, save_name)
            main_title = f'{tm.strftime("%Y-%m-%d %H:%M")} at {ht}m'

            if ht == 10:
                uctrl = np.squeeze(ds_ctrl['U10'])
                vctrl = np.squeeze(ds_ctrl['V10'])
            else:
                uctrl = np.squeeze(ds_ctrl.sel(height=ht)['U'])
                vctrl = np.squeeze(ds_ctrl.sel(height=ht)['V'])

            # calculate wind speed from u and v
            speed_ctrl = cf.wind_uv_to_spd(uctrl, vctrl)

            lon = speed_ctrl.XLONG.values
            lat = speed_ctrl.XLAT.values

            # calculate power from control file in kW
            # find the turbine location indices
            power_ctrl = np.array([])
            for i, row in turb_csv.iterrows():
                a = abs(lat - row.lat) + abs(lon - row.lon)
                i, j = np.unravel_index(a.argmin(), a.shape)
                power_calc = np.interp(speed_ctrl[i, j].values, power_curve['Wind Speed'], power_curve['Power'])
                power_ctrl = np.append(power_ctrl, power_calc)

            # calculate total windfarm power generated in GW
            cumulative_power_ctrl = np.round(np.sum(power_ctrl) / 1000000, 2)
            
            # add each hour for the 24-hour total
            cumulative_power_ctrl_total += cumulative_power_ctrl

            diff = speed_ctrl - speed_ctrl  # Since we're not plotting wind speed, set diff to zero

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
            cmap_power = plt.get_cmap('viridis')  # Choose an appropriate colormap for power
            levels_power = list(np.arange(power_min, power_max + power_interval, power_interval))
            color_label_power = 'Power (GW)'

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
            kwargs['cmap'] = cmap_power  # Use the power colormap
            kwargs['clab'] = color_label_power
            kwargs['levels'] = levels_power
            kwargs['extend'] = 'both'
            pf.plot_contourf(fig, ax1, lon, lat, power_ctrl, **kwargs)

            # plot difference (set to zeros since we're not plotting wind speed)
            kwargs['cmap'] = cmap_diff
            kwargs['clab'] = color_label_diff
            kwargs['levels'] = levels_diff
            kwargs['cbar_ticks'] = cticks
            pf.plot_contourf(fig, ax3, lon, lat, diff, **kwargs)

            # add WEA outline
            lease = glob.glob('/home/coolgroup/bpu/mapdata/shapefiles/BOEM-Renewable-Energy-Shapefiles-current/Wind_Lease_Outlines*.shp')[0]
            kwargs = dict()
            kwargs['edgecolor'] = 'magenta'
            pf.map_add_boem_outlines(ax1, lease, **kwargs)
            pf.map_add_boem_outlines(ax3, lease, **kwargs)

            # add turbine locations to control and diff
            ax1.scatter(turb_csv.lon, turb_csv.lat, s=.5, color='k', transform=ccrs.PlateCarree())
            ax3.scatter(turb_csv.lon, turb_csv.lat, s=.5, color='k', transform=ccrs.PlateCarree())

            ax1.set_xticklabels(ax1.get_xticklabels(), rotation=25, ha='center')
            ax2.set_xticklabels(ax2.get_xticklabels(), rotation=25, ha='center')
            ax3.set_xticklabels(ax3.get_xticklabels(), rotation=25, ha='center')

            ax1.set_title(f'Control (Power: {cumulative_power_ctrl} GW)', y=1.02)
            ax2.set_title(f'Wind Farm (Power: {cumulative_power_ctrl} GW)', y=1.02)
            ax3.set_title(f'Power Difference', y=1.02)
            plt.savefig(save_file, dpi=200)
            plt.close()

    # print statement for 24-hour power totals          
    print(f"Final Cumulative Control Power: {cumulative_power_ctrl_total} GW")


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

    arg_parser.add_argument('-power_min',
                            type=float,
                            default=0,
                            help='minimum colorbar limit for power plots')

    arg_parser.add_argument('-power_max',
                            type=float,
                            default=10,
                            help='maximum colorbar limit for power plots')

    arg_parser.add_argument('-power_interval',
                            type=float,
                            default=1,
                            help='interval for colorbar for power plots')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/windturbs/plots',
                            type=str,
                            help='Full directory path to save output plots.')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
