#!/usr/bin/env python

"""
Author: Lori Garzio on 4/25/2023
Last modified: Lori Garzio 9/13/2023
Grab U and V data at 10m and 160m from THREDDS for a user-defined time-range, calculate monthly averages and variance
and export as NetCDF.
"""

import argparse
import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
from collections import OrderedDict
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console


def main(args):
    start_str = args.start
    end_str = args.end
    domain = args.domain
    save_dir = args.save_dir
    loglevel = args.loglevel.upper()

    heights = [10, 160]

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    if domain == '3km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
    elif domain == '9km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_9km_processed/WRF_4.1_9km_Processed_Dataset_Best'
    elif domain == '1km_wf2km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_wf2km_processed/WRF_4.1_1km_with_Wind_Farm_Processed_Dataset_Best'
    elif domain == '1km_ctrl':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_ctrl_processed/WRF_4.1_1km_Control_Processed_Dataset_Best'
    else:
        raise ValueError('Invalid domain specified')

    save_dir = os.path.join(save_dir, domain)
    os.makedirs(save_dir, exist_ok=True)
    save_dir_logs = os.path.join(save_dir, 'logs')
    os.makedirs(save_dir_logs, exist_ok=True)

    # set up log file
    logfile = os.path.join(save_dir_logs, f'{domain}_monthly-avgs.log')
    logging = cf.setup_logger(f'logging_{domain}', loglevel, logfile)

    # check if a NetCDF file already exists
    save_file = os.path.join(save_dir, f'ruwrf_{domain}_monthly_averages.nc')
    if os.path.isfile(save_file):
        logging.warning(f'File for {domain} monthly averages already exists: {save_file}')
    else:
        ds = xr.open_dataset(mlink)
        ds = ds.sel(time=slice(start_date, end_date), height=160)
        #ds = ds.sel(time=test_date_range, height=160)
        ds = ds[["U", "V", "U10", "V10"]]

        months = np.unique(ds['time.month'])

        # create empty data dictionary to which data are appended
        data = {
            "coords": {
                "month": {"dims": "month", "data": np.array(months, dtype='int32')},
                "years": {"dims": "month", "data": np.array([], dtype='object')},
                "height": {"dims": "height", "data": np.array(heights, dtype='int32')},
                "XLAT": {"dims": ('south_north', 'west_east'), "data": ds.XLAT.values},
                "XLONG": {"dims": ('south_north', 'west_east'), "data": ds.XLONG.values}
            },
            "attrs": {
                "model_domain": domain,
            },
            "dims": {"month"},
            "data_vars": dict()
        }

        variables = ['ws_monthly_avg', 'ws_monthly_variance']
        for v in variables:
            data["data_vars"][v] = dict()
            data["data_vars"][v]["data"] = np.empty((len(months), len(heights), np.shape(ds.XLAT)[0], np.shape(ds.XLAT)[1]), dtype='float32')
            data["data_vars"][v]["data"][:] = np.nan
            data["data_vars"][v]["dims"] = ('month', 'height', 'south_north', 'west_east')
            data["data_vars"][v]["attrs"] = dict()

        # add attrs for lat/lon
        latlonvars = ['XLONG', 'XLAT']
        for llv in latlonvars:
            data["coords"][llv]["attrs"] = ds[llv].attrs

        logging.info(f'Calculating variables {variables} for {domain} monthly averages: {start_str} to {end_str}')

        for i, month in enumerate(months):
            # calculate monthly-averaged windspeed at 10m and 160m
            month_idx = ds['time.month'].values == month
            years = ' '.join(str(x) for x in np.unique(ds['time.year'].values[month_idx]))
            data['coords']['years']['data'] = np.append(data['coords']['years']['data'], years)
            for ih, height in enumerate(heights):
                if height == 10:
                    u_month = ds.U10[month_idx]
                    v_month = ds.V10[month_idx]
                else:
                    u_month = ds.U[month_idx]
                    v_month = ds.V[month_idx]
                ws_month = cf.wind_uv_to_spd(u_month, v_month)
                ws_monthly = ws_month.mean(dim='time')

                # append monthly-averaged data to data dictionary
                data["data_vars"]['ws_monthly_avg']['data'][i][ih] = ws_monthly
                data["data_vars"]['ws_monthly_avg']["attrs"]["units"] = 'm s-1'
                data["data_vars"]['ws_monthly_avg']["attrs"]["long_name"] = 'Monthly Averaged Wind Speed'
                data["data_vars"]['ws_monthly_avg']["attrs"]["comment"] = f'Average of wind speeds for each month for ' \
                                                                          f'specified years'

                # # calculate monthly variance (taking wind speed magnitude and direction into account)
                # u_month_variance = u_month.var(dim='time')
                # v_month_variance = v_month.var(dim='time')
                # rms = np.sqrt((u_month_variance + v_month_variance))
                #
                # # append variance to data dictionary
                # data["data_vars"]['rms_monthly']['data'][i] = rms
                # data["data_vars"]['rms_monthly']["attrs"]["comment"] = f'Measure of variance for wind speed for each month for ' \
                #                                                f'years: {years}, calculated as sqrt(u_variance + v_variance)'

                # calculate monthly variance (just wind speed magnitude)
                ws_monthly_variance = ws_month.var(dim='time')

                # append variance to data dictionary
                data["data_vars"]['ws_monthly_variance']['data'][i][ih] = ws_monthly_variance
                data["data_vars"]['ws_monthly_variance']["attrs"]["units"] = 'm2 s-2'
                data["data_vars"]['ws_monthly_variance']["attrs"]["comment"] = f'Variance of wind speed magnitude for each ' \
                                                                               f'month for specified years'

        outds = xr.Dataset.from_dict(data)

        # add created time to global attrs
        datetime_format = '%Y-%m-%dT%H:%M:%SZ'
        created = dt.datetime.utcnow().strftime(datetime_format)  # creation time Timestamp
        time_start = pd.to_datetime(np.nanmin(ds.time.values)).strftime(datetime_format)
        time_end = pd.to_datetime(np.nanmax(ds.time.values)).strftime(datetime_format)

        ga_comment = f'Summarized data from RU-WRF dataset. Average and variance of wind speeds for each month for ' \
                     f'multiple years spanning: {start_str} to {end_str}'

        global_attributes = OrderedDict([
            ('date_created', created),
            ('date_modified', created),
            ('time_coverage_start', time_start),
            ('time_coverage_end', time_end),
            ('comment', ga_comment),
            ('data_source', mlink),
            ('model_version', 'v4.1')
        ])

        global_attributes.update(outds.attrs)

        outds = outds.assign_attrs(global_attributes)

        # add attrs for height and years
        outds.height.attrs = ds.height.attrs
        outds.years.attrs = dict(comment='Years included in each monthly average')

        # Add compression to all variables
        encoding = {}
        for k in outds.data_vars:
            encoding[k] = {'zlib': True, 'complevel': 1}

        # save .nc file
        outds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4')
        logging.info(f'Finished calculating monthly averages for {domain} dataset, 10 and 160m: {start_str} to {end_str}')
        overall_start = pd.to_datetime(np.nanmin(ds.time.values)).strftime('%Y-%m-%dT%H:%M')
        overall_end = pd.to_datetime(np.nanmax(ds.time.values)).strftime('%Y-%m-%dT%H:%M')
        logging.info(f'Date range summarized in {save_file}: {overall_start} to {overall_end}')


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-s', '--start',
                            dest='start',
                            default='20190101',
                            type=str,
                            help='Start Date in format YYYYMMDD ')

    arg_parser.add_argument('-e', '--end',
                            dest='end',
                            default='20221231',
                            type=str,
                            help='End Date in format YYYYMMDD')

    arg_parser.add_argument('-d', '--domain',
                            dest='domain',
                            default='3km',
                            type=str,
                            choices=['3km', '9km', '1km_wf2km', '1km_ctrl'],
                            help='Domain: operational 3km or 9km WRF domain, 1km with simulated windfarm (1km_wf2km), ' \
                                 '1km without simulated windfarm (1km_ctrl)')

    arg_parser.add_argument('-save_dir',
                            default='/home/coolgroup/bpu/wrf/data/wrf_nc/averages',
                            type=str,
                            help='Full directory path to save output files.')

    arg_parser.add_argument('-loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
