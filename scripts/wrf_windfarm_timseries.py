# This code will grab data from the RUWRF Upwelling vs No upwelling runs and plot an hourly time series for the time range specified.
# We are going to grab the nearest turbine locations and calculte the power calculate the power difference

import argparse
import sys
import numpy as np
import pandas as pd
import os
import glob
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, DayLocator

def process_date(ymd,lease):

    #Next we need to grab the arguments
    heights = [160]

    ###### Type of data we are going to grab, for now we only need to focus on 2, full Lease upwelling control run and Full lease upwelling #######

    upwell_fname ='1km_wf2km_nyb'
    noupwell_fname = '1km_wf2km_nyb_modsst'


    power_curve = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_lw15mw_power.csv')

    ###### create file names that are going to be retreived on the server ########

    #wrf_turb_dir = '/home/coolgroup/ru-wrf/real-time/v4.1_parallel/processed_windturbs' # contains the main directory for processed wind turbine files
    wrf_turb_dir = '/home/coolgroup/ru-wrf/real-time/v4.1_parallel/processed_windturbs'


    # create file names that are going to be retreived on the server
    files_upwell = sorted(glob.glob(os.path.join(wrf_turb_dir,upwell_fname ,ymd, '*.nc'))) 

    # Initialize cumulative power totals
    cumulative_power_total = 0
    cumulative_power_ctrl_total = 0
    cumulative_power_diff_total = 0

    # Create empty lists to store cumulative power values and corresponding timestamps
    cumulative_power_values = []
    cumulative_power_ctrl_values = []
    cumulative_power_diff_values = []

    #empty time values
    time_values_list = []

    for fname in files_upwell:
        if fname.split('_')[-1] == 'H000.nc':
            continue
        f = fname.split('/')[-1]

        # find the corresponding control file
        fname_ctrl = os.path.join(wrf_turb_dir, noupwell_fname, ymd, f)

        ds = xr.open_dataset(fname)
        ds_ctrl = xr.open_dataset(fname_ctrl)

        # Assuming there's a time dimension in your datasets
        time_values = pd.to_datetime(ds.Time.values)

        for ht in heights:
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

            # Calculate wind speed from u and v
            speed = np.sqrt(u**2 + v**2)
            speed_ctrl = np.sqrt(uctrl**2 + vctrl**2)

            # Initialize power arrays
            power = np.zeros_like(ds.XLAT.values)
            power_ctrl = np.zeros_like(ds.XLAT.values)

            # Calculate power at each turbine location
            for i, row in lease.iterrows():
                a = np.abs(ds.XLAT.values - row.lat) + np.abs(ds.XLONG.values - row.lon)
                i, j = np.unravel_index(a.argmin(), a.shape)

                # Calculate power at each turbine location
                power_calc = np.interp(speed[i, j], power_curve['Wind Speed'], power_curve['Power'])
                power[i, j] = power_calc

                # Calculate power from control file at each turbine location
                power_ctrl_calc = np.interp(speed_ctrl[i, j], power_curve['Wind Speed'], power_curve['Power'])
                power_ctrl[i, j] = power_ctrl_calc

            # Calculate total wind farm power generated in GW
            cumulative_power = np.sum(power) / 1000000
            cumulative_power_ctrl = np.sum(power_ctrl) / 1000000
            cumulative_power_diff = cumulative_power - cumulative_power_ctrl

            # Add each hour for the total
            cumulative_power_total += cumulative_power
            cumulative_power_ctrl_total += cumulative_power_ctrl
            cumulative_power_diff_total += cumulative_power_diff

            # Append values to lists
            cumulative_power_values.append(cumulative_power)
            cumulative_power_ctrl_values.append(cumulative_power_ctrl)
            cumulative_power_diff_values.append(cumulative_power_diff)

            #append time values
            time_values_list.extend(time_values)

    return time_values_list,cumulative_power_values, cumulative_power_diff_values, cumulative_power_ctrl_values ,cumulative_power_total, cumulative_power_ctrl_total, cumulative_power_diff_total


def main(args):
    start = args.start
    end = args.end
    location = args.loc
    month = args.month

    if location == 'full':
       lease = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/turbine_locations_final_nyb.csv')
       label_plot = 'All Lease Areas'
    elif location == 'coastal':
       lease = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/turbine_locations_final.csv')
       label_plot = 'Coastal Lease Areas'
    else:
        raise ValueError(f'Invalid location provided')

    start_date = pd.to_datetime(start)
    end_date = pd.to_datetime(end)
    date_range = pd.date_range(start=start_date, end=end_date, freq='D')

    # Create empty lists to store cumulative power values and corresponding timestamps
    f_cumulative_power_values = []
    f_cumulative_power_ctrl_values = []
    f_cumulative_power_diff_values = []
    f_time_values_list = []

    f_cumulative_power_total = 0
    f_cumulative_power_ctrl_total = 0
    f_cumulative_power_diff_total = 0

    for date in date_range:
        ymd = date.strftime('%Y%m%d')
        time_vals, cum_power_vals, cum_power_diff_vals, cum_power_ctrl_vals, cum_power_total, cum_power_ctrl_total, cum_power_diff_total = process_date(ymd,lease)

        # Add each hour for the total
        f_cumulative_power_total += cum_power_total
        f_cumulative_power_ctrl_total += cum_power_ctrl_total
        f_cumulative_power_diff_total += cum_power_diff_total

        # Extend lists instead of appending them
        f_cumulative_power_values.extend(cum_power_vals)
        f_cumulative_power_ctrl_values.extend(cum_power_ctrl_vals)
        f_cumulative_power_diff_values.extend(cum_power_diff_vals)

        # Append time values
        f_time_values_list.extend(time_vals)


    # Plotting time series


    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(14, 8))

    ax1.plot(f_time_values_list, f_cumulative_power_values, label='Upwelling')
    ax1.plot(f_time_values_list, f_cumulative_power_ctrl_values, label='No Upwelling')
    ax1.set_ylabel('Cumulative Power (GW)')
    ax1.set_title(f'{label_plot} Windfarm Hourly Power Production ({month})')
    ax1.legend()  # Add legend for the first axis
    ax1.grid()

    ax2.plot(f_time_values_list, f_cumulative_power_diff_values, label='Cumulative Power Difference', color='orange')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Hourly Power Difference (GW)')
    ax2.legend()  # Add legend for the second axis
    ax2.grid()
    # Set x-axis ticks to display only the day of the month
    ax2.xaxis.set_major_locator(DayLocator())
    ax2.xaxis.set_major_formatter(DateFormatter('%m/%d/%y'))
    plt.xticks(rotation=45)


    filepath = "/www/web/rucool/windenergy/ru-wrf/images/upwelling_case/upwelling_turb_timeseries/"+f'{month}_timeseries.png'
    plt.savefig(filepath)



if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description='Calculate daily WRF SST input power statistics',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('start',
                            type=str,
                            help='Start date for the range of dates in the format YYYYmmdd (e.g., 20220601).')

    arg_parser.add_argument('end',
                            type=str,
                            help='End date for the range of dates in the format YYYYmmdd (e.g., 20220630).')

    arg_parser.add_argument('-loc',
                            type=str,
                            default='full',
                            choices=['full','coastal'],
                            help='turbine layout')
    
    arg_parser.add_argument('month',
                            type=str,
                            help='add the month of the date range specified')
    


    parsed_args = arg_parser.parse_args()
    main(parsed_args)