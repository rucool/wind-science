#!/usr/bin/env python

"""
Author: Lori Garzio on 2/2/2023
Last modified: Lori Garzio 2/3/2023
Grab RU-WRF Power, U and V data from THREDDS for a user-defined time-range and height at the simulated turbine
locations and export as NetCDF.
"""

import argparse
import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
from collections import OrderedDict
from functions.common import setup_logger
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console


def main(args):
    start_str = args.start
    end_str = args.end
    height = args.height
    turb_csv = args.turb_csv
    domain = args.domain
    save_dir = args.save_dir
    loglevel = args.loglevel.upper()

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    variables = ["U", "V"]

    if domain == '1km_wf2km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_wf2km_processed/WRF_4.1_1km_with_Wind_Farm_Processed_Dataset_Best'
        turbs_on_off = 'on'
        variables.append("POWER")
    elif domain == '1km_ctrl':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_1km_ctrl_processed/WRF_4.1_1km_Control_Processed_Dataset_Best'
        turbs_on_off = 'off'
    else:
        raise ValueError('Invalid domain specified')

    save_dir = os.path.join(save_dir, domain)
    os.makedirs(save_dir, exist_ok=True)
    save_dir_logs = os.path.join(save_dir, 'logs')
    os.makedirs(save_dir_logs, exist_ok=True)

    # set up log file
    logfile = os.path.join(save_dir_logs, f'{domain}_{height}.log')
    logging = setup_logger(f'logging_{domain}', loglevel, logfile)

    # check if a NetCDF file already exists
    save_file = os.path.join(save_dir, f'{domain}_{height}.nc')
    if os.path.isfile(save_file):
        logging.warning(f'Subset file for {domain} and height {height}m already exists: {save_file}')
        logging.warning('Please use the wrf_data_updater_windturbs.py code to add more data to the existing file')
    else:
        ds = xr.open_dataset(mlink)
        ds = ds.sel(time=slice(start_date, end_date), height=height)
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

        # create empty data dictionary to which data are appended
        data = {
            "coords": {
                "time": {"dims": "time", "data": np.array([], dtype='datetime64[ns]')},
                "points": {"dims": "points", "data": np.array([], dtype='int32')},
                "height": {"dims": (), "data": np.array(height, dtype='int32')},
                "XLAT": {"dims": "points", "data": np.array([], dtype='float32')},
                "XLONG": {"dims": "points", "data": np.array([], dtype='float32')}
            },
            "attrs": {
                "model_domain": domain,
                "wind_turbines": turbs_on_off
            },
            "dims": {"time"},
            "data_vars": dict()
        }

        var_attrs = ['units', 'long_name', 'standard_name', 'description']
        for v in variables:
            data["data_vars"][v] = dict()
            data["data_vars"][v]["data"] = np.empty((len(ds.time), len(points)), dtype='float32')
            data["data_vars"][v]["data"][:] = np.nan
            data["data_vars"][v]["dims"] = ("time", "points")
            data["data_vars"][v]["attrs"] = dict()
            for va in var_attrs:
                try:
                    data["data_vars"][v]["attrs"][va] = ds[v].attrs[va]
                except KeyError:
                    continue

        # add attrs for lat/lon
        latlonvars = ['XLONG', 'XLAT']
        for llv in latlonvars:
            data["coords"][llv]["attrs"] = dict()
            for va in var_attrs:
                try:
                    data["coords"][llv]["attrs"][va] = ds[llv].attrs[va]
                except KeyError:
                    continue

        logging.info(f'Downloading variables {variables} for {domain} {height}m: {start_str} to {end_str}')

        for i, tm in enumerate(ds.time.values):
            # grab the data from the turbine locations and append to the empty data dictionary
            new_ds = xr.concat([ds.sel(time=tm, south_north=xi, west_east=yi) for xi, yi in points], dim='points')
            for v in variables:
                data["data_vars"][v]["data"][i] = new_ds[v].values

            data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)

        data['coords']['points']['data'] = np.array(new_ds.points.values, dtype='int32')
        data['coords']['XLONG']['data'] = np.array(new_ds.XLONG.values, dtype='float32')
        data['coords']['XLAT']['data'] = np.array(new_ds.XLAT.values, dtype='float32')

        outds = xr.Dataset.from_dict(data)

        # add created time to global attrs
        datetime_format = '%Y-%m-%dT%H:%M:%SZ'
        created = dt.datetime.utcnow().strftime(datetime_format)  # creation time Timestamp
        time_start = pd.to_datetime(outds.time.values[0]).strftime(datetime_format)
        time_end = pd.to_datetime(outds.time.values[-1]).strftime(datetime_format)

        global_attributes = OrderedDict([
            ('date_created', created),
            ('date_modified', created),
            ('time_coverage_start', time_start),
            ('time_coverage_end', time_end),
            ('comment', 'Data subset from RU-WRF dataset at the simulated turbine locations'),
            ('data_source', mlink),
            ('model_version', 'v4.1')
        ])

        global_attributes.update(outds.attrs)

        outds = outds.assign_attrs(global_attributes)

        # Add compression to all variables
        encoding = {}
        for k in outds.data_vars:
            encoding[k] = {'zlib': True, 'complevel': 1}

        encoding['time'] = dict(units='seconds since 1970-01-01 00:00:00', calendar='gregorian', zlib=False,
                                _FillValue=False, dtype=np.double)

        # save .nc file
        outds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')
        logging.info(f'Finished downloading data for {domain} {height}m: {start_str} to {end_str}')
        overall_start = pd.to_datetime(np.nanmin(outds.time.values)).strftime('%Y-%m-%dT%H:%M')
        overall_end = pd.to_datetime(np.nanmax(outds.time.values)).strftime('%Y-%m-%dT%H:%M')
        logging.info(f'Data range available in {save_file}: {overall_start} to {overall_end}')


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-s', '--start',
                            dest='start',
                            default='20220301',
                            type=str,
                            help='Start Date in format YYYYMMDD ')

    arg_parser.add_argument('-e', '--end',
                            dest='end',
                            default='20230228',
                            type=str,
                            help='End Date in format YYYYMMDD')

    arg_parser.add_argument('-z', '--height',
                            dest='height',
                            default=160,
                            type=list,
                            help='Height in meters')

    arg_parser.add_argument('-tcsv', '--turbine_csv',
                            dest='turb_csv',
                            default='/home/wrfadmin/toolboxes/wind-science/files/turbine_locations_final.csv',
                            help='Optional csv file containing coordinates for simulated turbine locations. '
                                 'Default is /home/wrfadmin/toolboxes/wind-science/files/turbine_locations_final.csv')

    arg_parser.add_argument('-d', '--domain',
                            dest='domain',
                            default='1km_wf2km',
                            type=str,
                            choices=['1km_wf2km', '1km_ctrl'],
                            help='Domain: 1km with simulated windfarm (1km_wf2km), 1km without simulated windfarm (1km_ctrl)')

    arg_parser.add_argument('-save_dir',
                            default='/home/coolgroup/bpu/wrf/data/wrf_nc/',
                            type=str,
                            help='Full directory path to save output files.')

    arg_parser.add_argument('-loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
