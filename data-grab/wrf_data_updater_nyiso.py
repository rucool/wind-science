#!/usr/bin/env python

"""
Author: Lori Garzio on 2/22/2024
Last modified: Lori Garzio 2/22/2024
Add more data to U, V, windspeed and power in existing NetCDF files (generated using wrf_data_wrangler_nyiso.py) from
the RU-WRF THREDDS server through a user-defined end date and height=160m (hub height) at specified locations.
Power is calculated using the NYISO power curve
"""

import argparse
import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
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
    end_str = args.end
    lease = args.lease
    domain = args.domain
    file_dir = args.file_dir
    loglevel = args.loglevel.upper()

    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)
    h = 160  # height

    if domain == '3km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
    elif domain == '9km':
        mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_9km_processed/WRF_4.1_9km_Processed_Dataset_Best'
    else:
        raise ValueError('Invalid domain specified')

    file_dir = os.path.join(file_dir, domain)
    file_dir_logs = os.path.join(file_dir, 'logs')

    ds = xr.open_dataset(mlink)
    wrf_lat = ds['XLAT']
    wrf_lon = ds['XLONG']

    pc = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_nyiso15mw_power.csv')
    location_csv = '/home/wrfadmin/toolboxes/wind-science/files/lease_centroids.csv'
    loc_df = pd.read_csv(location_csv)
    loc_df['lease_code'] = loc_df['lease'].map(lambda x: x.split(' - ')[0].replace(' ', ''))
    df = loc_df[loc_df['lease_code'] == lease]
    if len(df) == 0:
        raise ValueError(
            f'Invalid lease code: {lease}. Please provide a valid lease code, found in ./files/lease_centroids.csv (i.e. "OCS-A0512")')

    for midlon, midlat, co, st, ll, lease_code in zip(df['long'], df['lat'], df['company'], df['state'], df['lease'], df['lease_code']):
        # set up log file
        logfile = os.path.join(file_dir_logs, f'{lease_code}_{domain}.log')
        logging = cf.setup_logger(f'logging_{lease_code}', loglevel, logfile)

        # find the .nc file to which new data are appended
        ncfilename = os.path.join(file_dir, f'{lease_code}_{domain}_nyiso_power.nc')

        logging.info(f'Attempting to update U and V in {ncfilename}')

        if not os.path.isfile(ncfilename):
            logging.warning(f'No file exists for {lease_code} {domain}')
            logging.warning('Please use the wrf_data_wrangler_nyiso.py code to generate the original file')
            continue

        # find the last time stamp in the existing dataset, and convert the dataset to a dictionary so the new
        # data can be added
        try:
            with xr.open_dataset(ncfilename) as ncfile:
                ncfile = ncfile.load()
        except OSError as e:
            logging.error('Error reading file {:s} ({:})'.format(ncfilename, e))
            continue
        except ValueError as e:
            logging.error('Error reading file {:s} ({:})'.format(ncfilename, e))
            continue

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
            continue

        ds_filter = ds.sel(time=slice(tstart, end_date))
        logging.info(f'Downloading additional data for {lease_code} {h}m: {tstart.strftime("%Y%m%d")} to {end_str}')

        # find the closest WRF coordinate
        d = haversine_dist(midlon, midlat, wrf_lon, wrf_lat)
        dmin = np.argwhere(d.data == np.min(d.data))
        di = dmin[0][0]
        dj = dmin[0][1]

        ds_filter = ds_filter.sel(south_north=di, west_east=dj)
        ds_filter = ds_filter.sel(height=h)

        # make sure you found the same coordinates and height
        if ds_filter.XLAT.values != ncfile.XLAT.values:
            logging.warning(f'XLAT does not match between two datasets: {ds_filter.XLAT.values} and {ncfile.XLAT.values}. Exiting')
            continue
        if ds_filter.XLONG.values != ncfile.XLONG.values:
            logging.warning(f'XLONG does not match between two datasets: {ds_filter.XLONG.values} and {ncfile.XLONG.values}. Exiting')
            continue
        if ds_filter.height.values != ncfile.height.values:
            logging.warning(f'Height does not match between two datasets: {ds_filter.height.values} and {ncfile.height.values}. Exiting')
            continue

        ds_filter_u = ds_filter.U.values
        ds_filter_v = ds_filter.V.values
        data['coords']['time']['data'] = np.append(data['coords']['time']['data'], ds_filter.time.values)
        data['data_vars']['U']['data'] = np.append(data['data_vars']['U']['data'], ds_filter_u)
        data['data_vars']['V']['data'] = np.append(data['data_vars']['V']['data'], ds_filter_v)

        # for tm in ds_filter.time.values:
        #     data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
        #     data['data_vars']['U']['data'] = np.append(data['data_vars']['U']['data'], ds_filter.sel(time=tm)['U'].values)
        #     data['data_vars']['V']['data'] = np.append(data['data_vars']['V']['data'], ds_filter.sel(time=tm)['V'].values)

        # calculate wind speed and power
        speed = cf.wind_uv_to_spd(ds_filter_u, ds_filter_v)
        power = np.interp(speed, pc['Wind Speed'], pc['Power'])  # power in kW

        data['data_vars']['windspeed']['data'] = np.append(data['data_vars']['windspeed']['data'], speed)
        data['data_vars']['power']['data'] = np.append(data['data_vars']['power']['data'], power)

        # add modified time and update time range in global attrs
        datetime_format = '%Y-%m-%dT%H:%M:%SZ'
        modified = dt.datetime.utcnow().strftime(datetime_format)  # modified time Timestamp
        time_end = pd.to_datetime(data['coords']['time']['data'][-1]).strftime(datetime_format)

        data['attrs']['date_modified'] = modified
        data['attrs']['time_coverage_end'] = time_end

        outds = xr.Dataset.from_dict(data)

        # Add compression to all variables
        encoding = {}
        for k in outds.data_vars:
            encoding[k] = {'zlib': True, 'complevel': 1}

        encoding['time'] = dict(units='seconds since 1970-01-01 00:00:00', calendar='gregorian', zlib=False,
                                _FillValue=False, dtype=np.double)

        # save .nc file
        outds.to_netcdf(ncfilename, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')
        logging.info(f'Finished updating windspeeds for {lease_code} {h}m: {tstart.strftime("%Y%m%d")} to {end_str}')
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

    arg_parser.add_argument('-file_dir',
                            default='/www/web/rucool/windenergy/ru-wrf/nc_file_subset/lease_centroids',
                            type=str,
                            help='Full directory path to existing .nc files.')

    arg_parser.add_argument('-loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
