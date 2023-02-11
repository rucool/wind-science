#!/usr/bin/env python

"""
Author: Lori Garzio on 2/8/2023
Last modified: Lori Garzio 2/8/2023
Update data in existing NetCDF files (generated using wrf_data_wrangler_windturbs.py) from the RU-WRF THREDDS server
through a user-defined time-range and height at the simulated turbine locations.
"""

import argparse
import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
from functions.common import setup_logger
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console


def main(args):
    end_str = args.end
    height = args.height
    turb_csv = args.turb_csv
    domain = args.domain
    file_dir = args.file_dir
    loglevel = args.loglevel.upper()

    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    if domain == '1km_wf2km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_wf2km_processed/WRF_4.1_1km_with_Wind_Farm_Processed_Dataset_Best'
    elif domain == '1km_ctrl':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_ctrl_processed/WRF_4.1_1km_Control_Processed_Dataset_Best'
    else:
        raise ValueError('Invalid domain specified')

    file_dir = os.path.join(file_dir, domain)
    file_dir_logs = os.path.join(file_dir, 'logs')

    # set up log file
    logfile = os.path.join(file_dir_logs, f'{domain}_{height}.log')
    logging = setup_logger(f'logging_{domain}', loglevel, logfile)

    # find the .nc file to which new data are appended
    ncfilename = os.path.join(file_dir, f'{domain}_{height}.nc')

    logging.info(f'Attempting to update data in in {ncfilename}')

    if not os.path.isfile(ncfilename):
        logging.warning(f'No file exists for {domain} {height}m')
        logging.warning('Please use the wrf_data_wrangler_windturbs.py code to generate the original file')
    else:
        # find the last time stamp in the existing dataset, and convert the dataset to a dictionary so the new
        # data can be added
        try:
            with xr.open_dataset(ncfilename) as ncfile:
                ncfile = ncfile.load()
        except OSError as e:
            logging.error('Error reading file {:s} ({:})'.format(ncfilename, e))
        except ValueError as e:
            logging.error('Error reading file {:s} ({:})'.format(ncfilename, e))

        variables = list(ncfile.data_vars)
        tf = pd.to_datetime(np.nanmax(ncfile.time.values))
        data = ncfile.to_dict()

        # make sure the timestamps have the correct dtype
        data['coords']['time']['data'] = np.array(data['coords']['time']['data'], dtype='datetime64[ns]')

        # slice the WRF dataset on the timestamps to append to the subset dataset
        tstart = tf + dt.timedelta(hours=1)
        if end_date <= tstart:
            logging.warning(f'Data file end time: {tf.strftime("%Y-%m-%dT%H:%M")}')
            logging.warning('Trying to append data that already exists in the file. '
                            f'Choose an end date greater than the one provided: {end_str}')
        else:
            logging.info(f'Downloading additional data for {domain} {height}m: {tstart.strftime("%Y%m%d")} to {end_str}')
            ds = xr.open_dataset(mlink)
            ds = ds.sel(time=slice(tstart, end_date), height=height)
            ds = ds[variables]
            wrf_lat = ds['XLAT']
            wrf_lon = ds['XLONG']

            # find the turbine location indices
            df = pd.read_csv(turb_csv)
            points = []
            for i, row in df.iterrows():
                a = abs(wrf_lat - row.lat) + abs(wrf_lon - row.lon)
                i, j = np.unravel_index(a.argmin(), a.shape)
                points.append((i, j))

            for i, tm in enumerate(ds.time.values):
                # grab the data from the turbine locations and append to the data dictionary
                new_ds = xr.concat([ds.sel(time=tm, south_north=xi, west_east=yi) for xi, yi in points], dim='points')
                data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
                for v in variables:
                    data["data_vars"][v]["data"] = np.concatenate((data["data_vars"][v]["data"], [new_ds[v].values]))

            # add modified time and update time range in global attrs
            datetime_format = '%Y-%m-%dT%H:%M:%SZ'
            modified = dt.datetime.utcnow().strftime(datetime_format)  # modified time Timestamp
            time_end = pd.to_datetime(data['coords']['time']['data'][-1]).strftime(datetime_format)

            data['attrs']['date_modified'] = modified
            data['attrs']['time_coverage_end'] = time_end

            # assign the correct data types
            data['coords']['points']['data'] = np.array(data['coords']['points']['data'], dtype='int32')
            data['coords']['XLAT']['data'] = np.array(data['coords']['XLAT']['data'], dtype='float32')
            data['coords']['XLONG']['data'] = np.array(data['coords']['XLONG']['data'], dtype='float32')

            for v in variables:
                data['data_vars'][v]['data'] = np.array(data['data_vars'][v]['data'], dtype='float32')

            outds = xr.Dataset.from_dict(data)

            # Add compression to all variables
            encoding = {}
            for k in outds.data_vars:
                encoding[k] = {'zlib': True, 'complevel': 1}

            encoding['time'] = dict(units='seconds since 1970-01-01 00:00:00', calendar='gregorian', zlib=False,
                                    _FillValue=False, dtype=np.double)

            # save .nc file
            outds.to_netcdf(ncfilename, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')
            logging.info(f'Finished updating windspeeds for {domain} {height}m: {tstart.strftime("%Y%m%d")} to {end_str}')
            overall_start = pd.to_datetime(np.nanmin(outds.time.values)).strftime('%Y-%m-%dT%H:%M')
            overall_end = pd.to_datetime(np.nanmax(outds.time.values)).strftime('%Y-%m-%dT%H:%M')
            logging.info(f'Data range available in {ncfilename}: {overall_start} to {overall_end}')


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-e', '--end',
                            dest='end',
                            default='20221130',
                            type=str,
                            help='End Date in format YYYYMMDD')

    arg_parser.add_argument('-z', '--height',
                            dest='height',
                            default=160,
                            type=list,
                            help='Height in meters')

    arg_parser.add_argument('-d', '--domain',
                            dest='domain',
                            default='1km_wf2km',
                            type=str,
                            choices=['1km_wf2km', '1km_ctrl'],
                            help='Domain: 1km with simulated windfarm (1km_wf2km), 1km without simulated windfarm (1km_ctrl)')

    arg_parser.add_argument('-tcsv', '--turbine_csv',
                            dest='turb_csv',
                            default='/home/wrfadmin/toolboxes/wind-science/files/turbine_locations_final.csv',
                            help='Optional csv file containing coordinates for simulated turbine locations. '
                                 'Default is /home/wrfadmin/toolboxes/wind-science/files/turbine_locations_final.csv')

    arg_parser.add_argument('-file_dir',
                            default='/home/coolgroup/bpu/wrf/data/wrf_nc/',
                            type=str,
                            help='Full directory path to existing .nc files.')

    arg_parser.add_argument('-loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
