#!/usr/bin/env python

"""
Author: Lori Garzio on 11/10/2022
Last modified: Lori Garzio 11/21/2022
Grab RU-WRF u, v data from thredds for a user-defined time-range and height at specified locations
"""

import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
import pickle
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


def wind_uv_to_spd(u, v):
    """
    Calculates the wind speed from the u and v wind components
    :param u: west/east direction (wind from the west is positive, from the east is negative)
    :param v: south/noth direction (wind from the south is positive, from the north is negative)
    :returns wspd: wind speed calculated from the u and v wind components
    """
    wspd = np.sqrt(np.square(u) + np.square(v))

    return wspd


def main(start_str, end_str, height, save_dir):
    start_date = dt.datetime.strptime(start_str, '%Y%m%d')
    end_date = dt.datetime.strptime(end_str, '%Y%m%d') + dt.timedelta(hours=23)

    mlink = 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
    ds = xr.open_dataset(mlink)
    ds = ds.sel(time=slice(start_date, end_date))
    wrf_lat = ds['XLAT']
    wrf_lon = ds['XLONG']

    location_csv = '/Users/garzio/Documents/repo/rucool/wind-science/capacity_factor_analysis/files/lease_centroids.csv'
    loc_df = pd.read_csv(location_csv)
    for i, row in loc_df.iterrows():
        if np.logical_or(row.state == 'New Jersey', row.state == 'NY/NJ'):
            data = dict(time=np.array([], dtype='datetime64[ns]'),
                        u=np.array([]),
                        v=np.array([]))

            data['company'] = row.company
            data['code'] = row.code
            data['state'] = row.state
            data['lease'] = row.lease
            lease_code = row.lease.split(' - ')[0].replace(' ', '')

            midlon = row['long']
            midlat = row['lat']

            # find the closest WRF coordinate
            d = haversine_dist(midlon, midlat, wrf_lon, wrf_lat)
            dmin = np.argwhere(d.data == np.min(d.data))
            di = dmin[0][0]
            dj = dmin[0][1]

            ds_point = ds.sel(south_north=di, west_east=dj)
            ds_point_height = ds_point.sel(height=height)

            data['wrf_lat'] = ds_point_height.XLAT.values
            data['wrf_lon'] = ds_point_height.XLONG.values
            data['height'] = height

            for tm in ds_point_height.time.values:
                print(tm)
                data['time'] = np.append(data['time'], tm)
                data['u'] = np.append(data['u'], ds_point_height.sel(time=tm)['U'].values)
                data['v'] = np.append(data['v'], ds_point_height.sel(time=tm)['V'].values)

            pickle_name = os.path.join(save_dir, f'{lease_code}-{start_str}_{end_str}.pickle')
            with open(pickle_name, 'wb') as handle:
                pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
            print(f'Finished downloading {row["code"]}\n')


if __name__ == '__main__':
    start = '20201201'
    end = '20211130'
    height = 160
    sDir = '/Users/garzio/Documents/rucool/bpu/wrf/capacity_factor_analysis'
    main(start, end, height, sDir)
