#!/usr/bin/env python

"""
Author: lnazzaro 10/2022
Last modified: lnazzaro 11/6/2022
Calculate wind shear and alpha based on difference between two heights.
"""

# https://www.engineeringtoolbox.com/wind-shear-d_1215.html

import argparse
import sys
import numpy as np
from netCDF4 import Dataset, num2date
import xarray as xr
import pandas as pd
from datetime import timedelta

def calculate_directional_shear(z1,u1,v1,z2,u2,v2):
    # sort of made up based on https://learningweather.psu.edu/node/93
    # still needs a lot of checking
    dz=np.subtract(z2,z1)
    du=np.subtract(u2-u1)
    dv=np.subtract(v2-v1)
    dspeed=np.sqrt(np.square(du)+np.square(dv))
    s1=np.sqrt(np.square(u1)+np.square(v1))
    s2=np.sqrt(np.square(u2)+np.square(v2))
    u1_unit=np.divide(u1,s1)
    v1_unit=np.divide(v1,s1)
    u2_unit=np.divide(u2,s2)
    v2_unit=np.divide(v2,s2)
    du_unit=np.subtract(u2_unit,u1_unit)
    dv_unit=np.subtract(v2_unit,v1_unit)
    ddir=np.arctan2(dv_unit,du_unit)*180/np.pi
    shear_speed=np.divide(dspeed,dz)
    shear_dir=np.divide(ddir,dz)
    return shear_speed,shear_dir

def calculate_speed_shear(z1,s1,z2,s2):
    # delta speed / delta z
    speed_shear=np.divide(np.subtract(s2,s1),np.subtract(z2,z1))
    return speed_shear

def calculate_alpha(z1,s1,z2,s2):
    # https://csl.noaa.gov/projects/lamar/windshearformula.html
    # u=(uref)*((z/zref)^α), zref usually 10m
    alpha = np.log(np.divide(s2,s1))/np.log(np.divide(z2,z1))
    return alpha

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

def read_wrf(m,t0,t1,lon,lat,h0,h1,fullrange):
    mlink = f'https://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_{m}_processed/WRF_4.1_{m}_Processed_Dataset_Best'
    wrf = xr.open_dataset(mlink)
    mlon = wrf['XLONG']
    mlat = wrf['XLAT']
    mz = wrf['height']
    mtime = wrf['time']
    if t0==t1:
        time = mtime[np.abs(mtime-t0)==np.min(np.abs(mtime-t0))]
    else:
        time = mtime[np.logical_and(mtime>=t0, mtime<=t1)]
    z = mz[np.logical_and(mz>=h0, mz<=h1)]
    if not fullrange:
        z = z[[0,-1]]
    if np.size(lon)==1:
        d = haversine_dist(lon, lat, mlon, mlat)
        dmin = np.argwhere(d.data==np.min(d.data))
        di = dmin[0][0]
        dj = dmin[0][1]
        lon = mlon[di, dj]
        lat = mlat[di, dj]
        u = np.nan*np.ones((len(time),len(z),1,1))
        v = np.nan*np.ones((len(time),len(z),1,1))
        for t in time:
            if len(z)<=2:
                for zi in z:
                    u[time==t,z==zi,0,0] = wrf['U'][mtime==t,mz==zi,di,dj]
                    v[time==t,z==zi,0,0] = wrf['V'][mtime==t,mz==zi,di,dj]
            else:
                u[time==t,:,0,0] = wrf['U'][mtime==t,np.logical_and(mz>=np.min(z),mz<=np.max(z)),di,dj]
                v[time==t,:,0,0] = wrf['V'][mtime==t,np.logical_and(mz>=np.min(z),mz<=np.max(z)),di,dj]
    else:
        d = np.logical_and(
            np.logical_and(mlon>=lon[0], mlon<=lon[1]),
            np.logical_and(mlat>=lat[0], mlat<=lat[1])
        )
        di = np.sum(d,1)>0
        dj = np.sum(d,0)>0
        lon = mlon[di,dj]
        lat = mlat[di,dj]
        u = np.nan*np.ones((len(time),len(z),np.sum(di.data),np.sum(dj.data)))
        v = np.nan*np.ones((len(time),len(z),np.sum(di.data),np.sum(dj.data)))
        for t in time:
            if len(z)<=2:
                for zi in z:
                    u[time==t,z==zi,:,:] = wrf['U'][mtime==t,mz==zi,di,dj]
                    v[time==t,z==zi,:,:] = wrf['V'][mtime==t,mz==zi,di,dj]
            else:
                u[time==t,:,:,:] = wrf['U'][mtime==t,np.logical_and(mz>=np.min(z),mz<=np.max(z)),di,dj]
                v[time==t,:,:,:] = wrf['V'][mtime==t,np.logical_and(mz>=np.min(z),mz<=np.max(z)),di,dj]
    wrf.close()
    return time, z, lon, lat, u, v

def main(args):
    source = args.data_source
    t1 = args.end_time
    if t1 == 'latest':
        t1=pd.Timestamp.now()
    else:
        t1=pd.to_datetime(t1)
    t0=args.start_time
    if t0.isnumeric() and int(t0)<93*24:
        t0=t1-timedelta(hours=t0)
    elif t0.isnumeric() and int(t0)<10000000:
        print('Please provide starting date less than 93 days before t1, or format yyyymmddTHHMM')
        return
    else:
        t0=pd.to_datetime(t0)
    h0=args.lower_layer
    h1=args.upper_layer
    if args.western_longitude:
        if not args.eastern_longitude or not args.southern_latitude or not args.northern_latitude:
            print('Please provide entire bounding box')
            return
        lon = [args.western_longitude, args.eastern_longitude]
        lat = [args.southern_latitude, args.northern_latitude]
    else:
        lon = args.longitude
        lat = args.latitude
    calc_shear = args.shear
    use_direction = args.direction
    calc_alpha = args. alpha
    layer_method = args.shear_layers

    if layer_method=='maxdiff':
        if use_direction:
            print('Capability not yet set up to account for direction in shear calculation for maximum difference option. Ignoring direction.')
            use_direction = False
    
    # get data
    if source[:3]=='wrf':
        if layer_method in ['iterative','maxdiff']:
            time, z, lon, lat, u, v = read_wrf(source[-3:],t0,t1,lon,lat,h0,h1,True)
        elif layer_method=='topbottom':
            time, z, lon, lat, u, v = read_wrf(source[-3:],t0,t1,lon,lat,h0,h1,False)
        else:
            print('Please choose method from [iterative, maxdiff, and topbottom]')
            return
        speed=np.sqrt(np.square(u)+np.square(v))
        if layer_method=='maxdiff':
            s1 = np.min(speed,axis=1)
            s2 = np.max(speed,axis=1)
        else:
            s1 = np.squeeze(speed[:,0,:,:])
            s2 = np.squeeze(speed[:,-1,:,:])
    
    # calculate shear and/or alpha
    shear = None
    direction_shear = None
    alpha = None

    if layer_method=='iterative':
        for i in range(len(z)-1):
            if calc_shear:
                shear = np.nan*np.ones(np.shape(speed))
                if use_direction:
                    direction_shear = np.nan*np.ones(np.shape(speed))
                    ss, sd = calculate_directional_shear(z[i],u[:,i,:,:],v[:,i,:,:],z[i+1],u[:,i+1,:,:],v[:,i+1,:,:])
                    shear[:,i,:,:] = ss
                    direction_shear[:,i,:,:] = sd
                else:
                    ss = calculate_speed_shear(z[i],speed[:,i,:,:],z[i+1],speed[:,i+1,:,:])
                    shear[:,i,:,:] = ss
            if calc_alpha:
                alpha = np.nan*np.ones(np.shape(speed))
                a = calculate_alpha(z[i],speed[:,i,:,:],z[i+1],speed[:,i+1,:,:])
                alpha[:,i,:,:] = a
    else:
        if calc_shear:
            if use_direction:
                shear, direction_shear = calculate_directional_shear(z[0],u[:,0,:,:],v[:,0,:,:],z[-1],u[:,-1,:,:],v[:,-1,:,:])
            else:
                shear = calculate_speed_shear(z[0],s1,z[-1],s2)
        if calc_alpha:
            alpha = calculate_alpha(z[0],s1,z[-1],s2)
    
    return time, z, lon, lat, shear, direction_shear, alpha

    

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # data source
    arg_parser.add_argument('-f', '--data_source',
                            help='data source (currently only WRF; may add buoys)',
                            choices=['wrf3km','wrf9km'],
                            default='wrf3km')
    # T0
    arg_parser.add_argument('-t0', '--start_time',
                            help='start time (yyyymmddTHHMM, or number of hours to go backwards from t1)',
                            default=0)
    # T1
    arg_parser.add_argument('-t1', '--end_time',
                            help='end time (yyyymmddTHHMM, or latest)',
                            default='latest')
    # H0
    arg_parser.add_argument('-h0', '--lower_layer',
                            help='lower height layer (m) used to calculate shear',
                            choices=[10,40,60,80,100,120,140,160,180,200,220,240,250,260,280,300,320],
                            default=10)
    # H1
    arg_parser.add_argument('-h1', '--upper_layer',
                            help='upper height layer (m) used to calculate shear',
                            choices=[10,40,60,80,100,120,140,160,180,200,220,240,250,260,280,300,320],
                            default=80)
    # location
    arg_parser.add_argument('-lon', '--longitude',
                            help='longitude (single data point)',
                            default=-74.08194)
    arg_parser.add_argument('-lat', '--latitude',
                            help='latitude (single data point)',
                            choices=['rt', 'delayed'],
                            default=39.2025)
    # min/max lon
    arg_parser.add_argument('-minlon', '--western_longitude',
                            help='minimum (westernmost) longitude for bounding box (empty if single point)',
                            default=None)
    arg_parser.add_argument('-maxlon', '--eastern_longitude',
                            help='maximum (easternmost) longitude for bounding box (empty if single point)',
                            default=None)
    # min/max lat
    arg_parser.add_argument('-minlat', '--southern_latitude',
                            help='minimum (southernmost) latitude for bounding box (empty if single point)',
                            default=None)
    arg_parser.add_argument('-maxlat', '--northern_latitude',
                            help='maximum (northernmost) longitude for bounding box (empty if single point)',
                            default=None)
    # shear/direction/alpha - which to calculate
    arg_parser.add_argument('-s', '--shear',
                            help='whether or not to calculate shear (T/F)',
                            choices=[True, False],
                            default=True)
    arg_parser.add_argument('-d', '--direction',
                            help='whether or not to account for direction change in shear calculation (T/F)',
                            choices=[True, False],
                            default=False)
    arg_parser.add_argument('-a', '--alpha',
                            help='whether or not to calculate alpha',
                            choices=[True, False],
                            default=False)
    # type of calculation across height
    arg_parser.add_argument('-t', '--shear_layers',
                            help='way to calculate shear across multiple layers (iterative: shear calculated between each neighboring height layer; maxdiff: shear based on maximum speed difference anywhere in given altitude range; topbottom: shear based on difference between upper layer and lower layer)',
                            choices=['iterative','maxdiff','topbottom'],
                            default='topbottom')

    parsed_args = arg_parser.parse_args()

    sys.exit(main(parsed_args))