import argparse
import sys
import numpy as np
import pandas as pd
import os
import glob
import xarray as xr

def main(args):
    start_date = pd.to_datetime(args.start_date)
    end_date = pd.to_datetime(args.end_date)
    date_range = pd.date_range(start=start_date, end=end_date, freq='D')

    cumulative_power_total = 0
    cumulative_power_ctrl_total = 0
    cumulative_power_diff_total = 0

    for date in date_range:
        ymd = date.strftime('%Y%m%d')
        args.ymd = ymd
        cumulative_power, cumulative_power_ctrl, cumulative_power_diff = process_date(args)
        cumulative_power_total += cumulative_power
        cumulative_power_ctrl_total += cumulative_power_ctrl
        cumulative_power_diff_total += cumulative_power_diff

    # print statement for total power statistics
    print(f"Total Cumulative Power: {cumulative_power_total} GW")
    print(f"Total Cumulative Control Power: {cumulative_power_ctrl_total} GW")
    print(f"Total Cumulative Power Difference: {cumulative_power_diff_total} GW")

def process_date(args):
    ymd = args.ymd
    expt = args.expt
    heights = args.heights

    turb_csv_dir = '/home/wrfadmin/toolboxes/wind-science/files'
    power_curve = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_lw15mw_power.csv')
    file_dir = '/home/coolgroup/ru-wrf/real-time/v4.1_parallel/processed_windturbs'

    control_dir = "1km_wf2km_nyb"  # Adjust this if needed

    # Load turbine CSV file
    if expt == '1km_wf2km':
        turb_csv = pd.read_csv(os.path.join(turb_csv_dir, 'turbine_locations_final.csv'))
    elif expt == '1km_wf2km_nyb' or expt == '1km_wf2km_nyb_modsst':
        turb_csv = pd.read_csv(os.path.join(turb_csv_dir, 'turbine_locations_final_nyb.csv'))
    else:
        raise ValueError(f'Invalid experiment provided: {expt}')

    files = sorted(glob.glob(os.path.join(file_dir, expt, ymd, '*.nc')))

    cumulative_power_total = 0
    cumulative_power_ctrl_total = 0
    cumulative_power_diff_total = 0

    for fname in files:
        if fname.split('_')[-1] == 'H000.nc':
            continue
        f = fname.split('/')[-1]

        # find the corresponding control file
        fname_ctrl = os.path.join(file_dir, control_dir, ymd, f)

        ds = xr.open_dataset(fname)
        ds_ctrl = xr.open_dataset(fname_ctrl)

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

            # calculate wind speed from u and v
            speed = np.sqrt(u**2 + v**2)
            speed_ctrl = np.sqrt(uctrl**2 + vctrl**2)

            # get power from wind farm file in kW
            power = ds.POWER.values / 1000

            # calculate power from control file in kW
            # find the turbine location indices
            power_ctrl = np.array([])
            for i, row in turb_csv.iterrows():
                a = np.abs(ds.XLAT.values - row.lat) + np.abs(ds.XLONG.values - row.lon)
                i, j = np.unravel_index(a.argmin(), a.shape)
                power_calc = np.interp(speed_ctrl[i, j], power_curve['Wind Speed'], power_curve['Power'])
                power_ctrl = np.append(power_ctrl, power_calc)

            # calculate total wind farm power generated in GW
            cumulative_power = np.sum(power) / 1000000
            cumulative_power_ctrl = np.sum(power_ctrl) / 1000000
            cumulative_power_diff = cumulative_power - cumulative_power_ctrl

            # add each hour for the total
            cumulative_power_total += cumulative_power
            cumulative_power_ctrl_total += cumulative_power_ctrl
            cumulative_power_diff_total += cumulative_power_diff

    return cumulative_power_total, cumulative_power_ctrl_total, cumulative_power_diff_total

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description='Calculate daily WRF SST input power statistics',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('start_date',
                            type=str,
                            help='Start date for the range of dates in the format YYYYmmdd (e.g., 20220601).')

    arg_parser.add_argument('end_date',
                            type=str,
                            help='End date for the range of dates in the format YYYYmmdd (e.g., 20220630).')

    arg_parser.add_argument('-expt',
                            type=str,
                            default='1km_wf2km',
                            choices=['1km_wf2km', '1km_wf2km_nyb', '1km_wf2km_nyb_modsst'],
                            help='Experiment version to calculate power statistics against the control.')

    arg_parser.add_argument('-heights',
                            type=list,
                            default=[160],
                            choices=[[160], [10, 160]],
                            help='List of heights to calculate power statistics.')

    parsed_args = arg_parser.parse_args()
    main(parsed_args)
