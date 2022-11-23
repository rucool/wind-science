#!/usr/bin/env python

"""
Author: Lori Garzio on 11/10/2022
Last modified: Lori Garzio 11/22/2022
Grab RU-WRF U, V data from THREDDS for a user-defined time-range and height at specified locations and export as NetCDF
"""

import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
from collections import OrderedDict
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


def main(start_str, end_str, height, lease, state, save_dir):
    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=slice(start_date, end_date))
    wrf_lat = ds['XLAT']
    wrf_lon = ds['XLONG']

    location_csv = '/Users/garzio/Documents/repo/rucool/wind-science/capacity_factor_analysis/files/lease_centroids.csv'
    loc_df = pd.read_csv(location_csv)
    loc_df['lease_code'] = loc_df['lease'].map(lambda x: x.split(' - ')[0].replace(' ', ''))
    if lease:
        df = loc_df[loc_df['lease_code'] == lease]
    elif state:
        df = loc_df[loc_df['state'] == state]

    for midlon, midlat, co, st, ll, lease_code in zip(df['long'], df['lat'], df['company'], df['state'], df['lease'], df['lease_code']):
        # check if a NetCDF file already exists for the lease/height
        save_file = os.path.join(save_dir, f'{lease_code}_{height}.nc')
        if os.path.isfile(save_file):
            print(f'Subset file for this lease area and height already exists: {save_file}')
            print('Please use the wrf_data_updater.py code to add more data to the existing file')
            continue

        print(f'Downloading data for {lease_code} {height}m')

        # find the closest WRF coordinate
        d = haversine_dist(midlon, midlat, wrf_lon, wrf_lat)
        dmin = np.argwhere(d.data == np.min(d.data))
        di = dmin[0][0]
        dj = dmin[0][1]

        ds_filter = ds.sel(south_north=di, west_east=dj)
        ds_filter = ds_filter.sel(height=height)

        wrf_xlat = ds_filter.XLAT.values
        wrf_xlon = ds_filter.XLONG.values

        # create empty dictionary to which data are added
        data = {
            "coords": {
                "time": {"dims": "time", "data": np.array([], dtype='datetime64[ns]')},
                "height": {"dims": (), "data": np.array(height, dtype='int32')},
                "XLAT": {"dims": (), "data": np.array(wrf_xlat, dtype='float32')},
                "XLONG": {"dims": (), "data": np.array(wrf_xlon, dtype='float32')}
            },
            "attrs": {
                "point_lat": midlat,
                "point_lon": midlon,
                "company": co,
                "state": st,
                "lease": ll,
                "lease_code": lease_code
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
                                   f"XLONG {str(np.round(wrf_xlon, 6))}, height {height}m"
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
                                   f"XLONG {str(np.round(wrf_xlon, 6))}, height {height}m"
                    }
                },
            },
        }

        for tm in ds_filter.time.values:
            data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
            data['data_vars']['U']['data'] = np.append(data['data_vars']['U']['data'], ds_filter.sel(time=tm)['U'].values)
            data['data_vars']['V']['data'] = np.append(data['data_vars']['V']['data'], ds_filter.sel(time=tm)['V'].values)

        outds = xr.Dataset.from_dict(data)

        # add created time to global attrs
        datetime_format = '%Y-%m-%dT%H:%M:%SZ'
        created = dt.datetime.utcnow().strftime(datetime_format)  # creation time Timestamp
        time_start = pd.to_datetime(outds.time.values[0]).strftime(datetime_format)
        time_end = pd.to_datetime(outds.time.values[-1]).strftime(datetime_format)

        global_attributes = OrderedDict([
            ('date_created', created),
            ('time_coverage_start', time_start),
            ('time_coverage_end', time_end),
            ('comment', 'U and V subset from RU-WRF dataset at the specified location and height'),
            ('data_source', mlink)
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
        print(f'Finished downloading {lease_code} {height}m')


if __name__ == '__main__':
    start = '20190101'
    end = '20190101'
    height = 160
    lease = None  # 'OCS-A0498'
    state = 'New Jersey'  # 'New Jersey', 'NY/NJ'
    sDir = '/Users/garzio/Documents/repo/rucool/wind-science/capacity_factor_analysis/files'
    main(start, end, height, lease, state, sDir)
