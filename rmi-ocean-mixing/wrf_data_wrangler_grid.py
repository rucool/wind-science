#!/usr/bin/env python

"""
Author: Lori Garzio on 4/25/2023
Last modified: Lori Garzio 9/25/2024
Grab WRF data (entire domain) from THREDDS for a user-defined time-range and export as NetCDF.
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
    heights = args.heights
    save_dir = args.save_dir
    loglevel = args.loglevel.upper()

    file_names = ['T2', 'TEMP', 'ws', 'UST', 'Q2', 'TKE_PBL', 'HFX']
    subset_vars = [['T2'], ['TEMP'], ['U', 'V', 'U10', 'V10'], ['UST'], ['Q2'], ['TKE_PBL'], ['HFX']]

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
    elif domain == '1km_wf2km_nyb':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/case_studies/wrf_4_1_1km_wf2km_upwellingctrl_processed/WRF_4.1_1km_with_Wind_Farm_NY_Bight_Control_Processed_Dataset_Best'
    else:
        raise ValueError('Invalid domain specified')

    save_dir = os.path.join(save_dir, domain)
    os.makedirs(save_dir, exist_ok=True)

    # set up log file
    logfile = os.path.join(save_dir, f'{domain}_{start_str}_{end_str}_subsetting.log')
    logging = cf.setup_logger(f'logging_{domain}', loglevel, logfile)

    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=slice(start_date, end_date), height=heights)
    for i, sv in enumerate(subset_vars):
        # subset data for each variable group
        save_file = os.path.join(save_dir, f'ruwrf_{domain}_{file_names[i]}_{start_str}_{end_str}.nc')
        logging.info(f'Subsetting WRF THREDDS {domain} variables: {sv} {start_str} to {end_str}')
        ds2 = ds[sv]

        # add created time to global attrs
        datetime_format = '%Y-%m-%dT%H:%M:%SZ'
        created = dt.datetime.utcnow().strftime(datetime_format)  # creation time Timestamp
        time_start = pd.to_datetime(np.nanmin(ds2.time.values)).strftime(datetime_format)
        time_end = pd.to_datetime(np.nanmax(ds2.time.values)).strftime(datetime_format)

        ga_comment = f'Subset model output from RU-WRF {domain} dataset from {start_str} to {end_str}'

        global_attributes = OrderedDict([
            ('date_created', created),
            ('date_modified', created),
            ('time_coverage_start', time_start),
            ('time_coverage_end', time_end),
            ('comment', ga_comment),
            ('data_source', mlink),
            ('model_version', 'RU-WRF v4.1'),
            ('creator_url', ds.creator_url),
            ('creator_name', ds.creator_name),
            ('creator_type', ds.creator_type),
            ('creator_email', ds.creator_email),
            ('creator_institution', ds.creator_institution),
            ('institution', ds.institution),
            ('acknowledgement', ds.acknowledgement),
            ('project', ds.project),
            ('geospatial_lat_min', ds.geospatial_lat_min),
            ('geospatial_lat_max', ds.geospatial_lat_max),
            ('geospatial_lon_min', ds.geospatial_lon_min),
            ('geospatial_lon_max', ds.geospatial_lon_max),
            ('geospatial_lat_units', ds.geospatial_lat_units),
            ('geospatial_lon_units', ds.geospatial_lon_units),
            ('contributor_name', ds.contributor_name),
            ('contributor_role', ds.contributor_role),
            ('references', 'https://rucool.marine.rutgers.edu/research/offshore-wind/'),
            ('MAP_PROJ_CHAR', ds.MAP_PROJ_CHAR),
            ('CEN_LAT', ds.CEN_LAT),
            ('CEN_LON', ds.CEN_LON),
            ('TRUELAT1', ds.TRUELAT1),
            ('TRUELAT2', ds.TRUELAT2),
            ('MOAD_CEN_LAT', ds.MOAD_CEN_LAT),
            ('STAND_LON', ds.STAND_LON),
            ('DX', ds.DX),
            ('DY', ds.DY)
        ])

        global_attributes.update(ds2.attrs)

        ds2 = ds2.assign_attrs(global_attributes)

        # Add compression to all variables
        encoding = {}
        for k in ds2.data_vars:
            # print(k)
            # edict = ds2[k].encoding
            # edict.update({'zlib': True, 'complevel': 1})
            # encoding[k] = edict
            encoding[k] = {'zlib': True, 'complevel': 1}

        # save .nc file
        ds2.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4')
        logging.info(f'Finished subsetting {domain} variables: {sv} {start_str} to {end_str}')
        overall_start = pd.to_datetime(np.nanmin(ds2.time.values)).strftime('%Y-%m-%dT%H:%M')
        overall_end = pd.to_datetime(np.nanmax(ds2.time.values)).strftime('%Y-%m-%dT%H:%M')
        logging.info(f'Date range in {save_file}: {overall_start} to {overall_end}')

    logging.info('Process finished')


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-s', '--start',
                            dest='start',
                            default='20220601',
                            type=str,
                            help='Start Date in format YYYYMMDD ')

    arg_parser.add_argument('-e', '--end',
                            dest='end',
                            default='20220630',
                            type=str,
                            help='End Date in format YYYYMMDD')

    arg_parser.add_argument('-d', '--domain',
                            dest='domain',
                            default='1km_wf2km_nyb',
                            type=str,
                            choices=['3km', '9km', '1km_wf2km', '1km_ctrl', '1km_wf2km_nyb'],
                            help='Domain: operational 3km or 9km WRF domain, 1km with simulated windfarm (1km_wf2km), ' \
                                 '1km without simulated windfarm (1km_ctrl)')

    arg_parser.add_argument('--heights',
                            dest='heights',
                            default=[120, 160, 200, 300],
                            type=list,
                            help='List of heights at which to extract data. Note the surface is automatically ' \
                                 'included.')

    arg_parser.add_argument('-save_dir',
                            type=str,
                            help='Full directory path to save output files.')

    arg_parser.add_argument('-loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
