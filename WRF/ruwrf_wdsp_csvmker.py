import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
import argparse
from functions import common as cf


# Sample command line: python ruwrf_wdsp_csvmaker.py 2022 1 1 2022 1 3 160


# Set up the argument parser
parser = argparse.ArgumentParser(description='Process WRF data.')
parser.add_argument('start_year', type=int, help='Start year (e.g., 2022)')
parser.add_argument('start_month', type=int, help='Start month (e.g., 1)')
parser.add_argument('start_day', type=int, help='Start day (e.g., 1)')
parser.add_argument('end_year', type=int, help='End year (e.g., 2022)')
parser.add_argument('end_month', type=int, help='End month (e.g., 1)')
parser.add_argument('end_day', type=int, help='End day (e.g., 3)')
parser.add_argument('height', type=int, help='Height (e.g., 160)')
args = parser.parse_args()

# Use the arguments to set the start and end dates and height
start_date = dt.datetime(args.start_year, args.start_month, args.start_day, 0, 0)  
end_date = dt.datetime(args.end_year, args.end_month, args.end_day, 23, 0)
height = args.height
point = [39.33202, -73.94720]

wlink = 'http://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
ds = xr.open_dataset(wlink)
ds = ds.sel(time=slice(start_date, end_date))

# Computes all distances between the wrf gridpoints for both lat and lon
a = abs(ds['XLAT'] - point[0]) + abs(ds['XLONG'] - point[1])
# calulcates the minimum values in a, which corresponds to the closest grid point and then converts the flat index to a 2d index
i, j = np.unravel_index(a.argmin(), a.shape)


########### Extract variables at the point ########

# get u and v component, calculate windspeed and estimated wind power
u = ds.sel(height = height)['U'][:, i, j].values
v = ds.sel(height = height)['V'][:, i, j].values
lat = ds['XLAT'][i, j].values.item()  # Get the scalar value
lon = ds['XLONG'][i, j].values.item()
times = ds['time'].values

# Calcualte wind speed and Direction
wind_speed = ws = cf.wind_uv_to_spd(u,v)
wind_direction = np.arctan2(-u, -v) * (180/np.pi)
wind_direction = (wind_direction + 360) % 360


# Assume 'times' is your numpy array of datetime64 values
dates = pd.to_datetime(times)

# Localize the datetime to UTC
dates_utc = dates.tz_localize('UTC')

# Create a DataFrame
df = pd.DataFrame({
    'Latitude': lat,
    'Longitude': lon,
    'Height (m)': height,
    'Date (UTC)': times,
    'U_velocity (m/s)': u,
    'V_velocity (m/s)': v,
    'Wind_Speed (m/s)': wind_speed,
    'Wind_Dir (degrees)': wind_direction
})

# Create the CSV filename with year and height
#csv_filename = f'/users/jameskim/Documents/rucool/Data/RUWRF_wind_data_{start_date.year}_{height}m.csv'
csv_filename = f'/home/wrfadmin/toolboxes/wind-science/files/RUWRF_wind_data_{start_date.year}_{height}m.csv'
df.to_csv(csv_filename, index=False)

print(f'Data saved to {csv_filename}')