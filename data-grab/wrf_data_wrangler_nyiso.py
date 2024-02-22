#!/usr/bin/env python

"""
Author: Lori Garzio on 2/19/2024
Last modified: Lori Garzio 2/21/2024
Grab RU-WRF U, V data from THREDDS for a user-defined time-range at height=160m (hub height) at specified locations.
Calculate windspeed and power using the NYISO power curve, and export as NetCDF.
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


# define equation to get distance (km) from one lon/lat to another (or an entire set)
def haversine_dist(blon,blat,slon,slat):
    # blon: longitude of single point
    # blat: latitude of single point
    # slon: longitude(s) of grid
    # slat: latitude(s) of grid
    R = 6373.0
    blon=blon*np.pi/180
    blat=blat*np.pi/180
    slon=slon*np.pi/180
    slat=slat*np.pi/180
    dlon=slon-blon
    dlat=slat-blat
    a=np.sin(dlat/2)**2+np.cos(blat)*np.cos(slat)*np.sin(dlon/2)**2
    c=2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
    distance=R*c
    return distance


def main(args):
    start_str = args.start
    end_str = args.end
    lease = args.lease
    domain = args.domain
    save_dir = args.save_dir
    loglevel = args.loglevel.upper()

    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)
    h = 160  # height

    if domain == '3km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
    elif domain == '9km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_9km_processed/WRF_4.1_9km_Processed_Dataset_Best'
    else:
        raise ValueError('Invalid domain specified')

    save_dir = os.path.join(save_dir, domain)
    os.makedirs(save_dir, exist_ok=True)
    save_dir_logs = os.path.join(save_dir, 'logs')
    os.makedirs(save_dir_logs, exist_ok=True)

    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=slice(start_date, end_date))
    wrf_lat = ds['XLAT']
    wrf_lon = ds['XLONG']

    pc = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_nyiso15mw_power.csv')
    location_csv = '/home/wrfadmin/toolboxes/wind-science/files/lease_centroids.csv'
    loc_df = pd.read_csv(location_csv)
    loc_df['lease_code'] = loc_df['lease'].map(lambda x: x.split(' - ')[0].replace(' ', ''))
    df = loc_df[loc_df['lease_code'] == lease]
    if len(df) == 0:
        raise ValueError(f'Invalid lease code: {lease}. Please provide a valid lease code, found in ./files/lease_centroids.csv (i.e. "OCS-A0512")')

    for midlon, midlat, co, st, ll, lease_code in zip(df['long'], df['lat'], df['company'], df['state'], df['lease'], df['lease_code']):
        # set up log file
        logfile = os.path.join(save_dir_logs, f'{lease_code}_{domain}.log')
        logging = cf.setup_logger(f'logging_{lease_code}', loglevel, logfile)

        logging.info(f'Attempting to download WRF data for state: {st}; company: {co}; lease {lease_code}')

        # check if a NetCDF file already exists for the lease/height
        save_file = os.path.join(save_dir, f'{lease_code}_{domain}_nyiso_power.nc')
        if os.path.isfile(save_file):
            logging.warning(f'Subset file for {lease_code} already exists: {save_file}')
            logging.warning('Please use the wrf_data_updater_nyiso.py code to add more data to the existing file')
            continue

        logging.info(f'Downloading WRF data for state: {st}; company: {co}; lease {lease_code}: {start_str} to {end_str}')

        # find the closest WRF coordinate
        d = haversine_dist(midlon, midlat, wrf_lon, wrf_lat)
        dmin = np.argwhere(d.data == np.min(d.data))
        di = dmin[0][0]
        dj = dmin[0][1]

        ds_filter = ds.sel(south_north=di, west_east=dj)
        ds_filter = ds_filter.sel(height=h)

        wrf_xlat = ds_filter.XLAT.values
        wrf_xlon = ds_filter.XLONG.values

        # create empty dictionary to which data are added
        data = {
            "coords": {
                "time": {"dims": "time", "data": np.array([], dtype='datetime64[ns]')},
                "height": {"dims": (), "data": np.array(h, dtype='int32')},
                "XLAT": {"dims": (), "data": np.array(wrf_xlat, dtype='float32')},
                "XLONG": {"dims": (), "data": np.array(wrf_xlon, dtype='float32')}
            },
            "attrs": {
                "point_lat": midlat,
                "point_lon": midlon,
                "company": co,
                "state": st,
                "lease": ll,
                "lease_code": lease_code,
                "model_domain": domain
            },
            "dims": "time",
            "data_vars": {
                "U": {
                    "dims": "time",
                    "data": np.array([], dtype='float32'),
                    "attrs": {
                        "units": "m s-1",
                        "long_name": "Eastward Wind Component",
                        "standard_name": "eastward_wind",
                        "description": "earth rotated u",
                        "comment": f"subset from RU-WRF at XLAT {str(np.round(wrf_xlat, 6))}, "
                                   f"XLONG {str(np.round(wrf_xlon, 6))}, height {h}m"
                    }
                },
                "V": {
                    "dims": "time",
                    "data": np.array([], dtype='float32'),
                    "attrs": {
                        "units": "m s-1",
                        "long_name": "Northward Wind Component",
                        "standard_name": "northward_wind",
                        "description": "earth rotated v",
                        "comment": f"subset from RU-WRF at XLAT {str(np.round(wrf_xlat, 6))}, "
                                   f"XLONG {str(np.round(wrf_xlon, 6))}, height {h}m"
                    }
                },
                "windspeed": {
                    "dims": "time",
                    "data": np.array([], dtype='float32'),
                    "attrs": {
                        "units": "m s-1",
                        "long_name": "Wind Speed",
                        "comment": f"Wind speed calculated from U and V subset from RU-WRF at XLAT "
                                   f"{str(np.round(wrf_xlat, 6))}, XLONG {str(np.round(wrf_xlon, 6))}, height {h}m"
                    }
                },
                "power": {
                    "dims": "time",
                    "data": np.array([], dtype='float32'),
                    "attrs": {
                        "units": "kW",
                        "long_name": "Estimated Wind Power",
                        "comment": f"Estimated wind power, calculated from wind speeds at at height of {h}m using"
                                   f"the NYSIO power curve."
                    }
                },
            },
        }

        data['coords']['time']['data'] = np.append(data['coords']['time']['data'], ds_filter.time.values)
        data['data_vars']['U']['data'] = np.append(data['data_vars']['U']['data'], ds_filter.U.values)
        data['data_vars']['V']['data'] = np.append(data['data_vars']['V']['data'], ds_filter.V.values)

        # for tm in ds_filter.time.values:
        #     data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
        #     data['data_vars']['U']['data'] = np.append(data['data_vars']['U']['data'], ds_filter.sel(time=tm)['U'].values)
        #     data['data_vars']['V']['data'] = np.append(data['data_vars']['V']['data'], ds_filter.sel(time=tm)['V'].values)

        # calculate wind speed and power
        speed = cf.wind_uv_to_spd(data['data_vars']['U']['data'], data['data_vars']['V']['data'])
        power = np.interp(speed, pc['Wind Speed'], pc['Power'])  # power in kW

        data['data_vars']['windspeed']['data'] = speed
        data['data_vars']['power']['data'] = power

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
            ('comment', 'U and V subset from RU-WRF dataset at the specified location and height. Windspeed is calculated'
                        'from U and V, and power is estimated from wind speed using the NYISO power curve.'),
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
        logging.info(f'Finished downloading windspeed data for {lease_code} {h}m {domain}: {start_str} to {end_str}')
        overall_start = pd.to_datetime(np.nanmin(outds.time.values)).strftime('%Y-%m-%dT%H:%M')
        overall_end = pd.to_datetime(np.nanmax(outds.time.values)).strftime('%Y-%m-%dT%H:%M')
        logging.info(f'Data range available in {save_file}: {overall_start} to {overall_end}')


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-s', '--start',
                            dest='start',
                            default='20190101',
                            type=str,
                            help='Start Date in format YYYYMMDD ')

    arg_parser.add_argument('-e', '--end',
                            dest='end',
                            default='20221031',
                            type=str,
                            help='End Date in format YYYYMMDD')

    arg_parser.add_argument('-l', '--lease',
                            dest='lease',
                            default='OCS-A0517',
                            help='Lease area, valid OCS lease codes can be found in '
                                 './files/lease_centroids.csv. Example: "OCS-A0512"')

    arg_parser.add_argument('-d', '--domain',
                            dest='domain',
                            default='3km',
                            type=str,
                            choices=['3km', '9km'],
                            help='Domain: operational 3km, operational 9km')

    arg_parser.add_argument('-save_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/nc_file_subset/lease_centroids',
                            type=str,
                            help='Full directory path to save output files.')

    arg_parser.add_argument('-loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
