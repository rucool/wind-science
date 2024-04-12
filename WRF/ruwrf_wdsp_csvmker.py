import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
import argparse
import os
from functions import common as cf


# Sample command line: python ruwrf_wdsp_csvmaker.py 2022 1 1 2022 1 3 160 "OCS-A 0499"


# Set up the argument parser
parser = argparse.ArgumentParser(description='Process WRF data.')
parser.add_argument('start_year', type=int, help='Start year (e.g., 2022)')
parser.add_argument('start_month', type=int, help='Start month (e.g., 1)')
parser.add_argument('start_day', type=int, help='Start day (e.g., 1)')
parser.add_argument('end_year', type=int, help='End year (e.g., 2022)')
parser.add_argument('end_month', type=int, help='End month (e.g., 1)')
parser.add_argument('end_day', type=int, help='End day (e.g., 3)')
parser.add_argument('height', type=int, help='Height (e.g., 160)')
parser.add_argument('location',type = str,help = 'Lease site: ( e.g., OCS-A 0499')
args = parser.parse_args()



def extract_coordinates(lease, lease_site):
    # Filter rows based on the specified site code
    selected_rows = lease[lease['lease'].str.contains(lease_site, case=False, regex=False)]

    # Extract latitude and longitude
    if not selected_rows.empty:
        latitude = selected_rows['lat'].iloc[0]
        longitude = selected_rows['long'].iloc[0]
        point = [latitude,longitude]
        return point
    else:
        print(f"No data found for the site code: {lease_site}")
        return None  # You can return None or any other value to indicate no data found



# Use the arguments to set the start and end dates and height
start_date = dt.datetime(args.start_year, args.start_month, args.start_day, 0, 0)  
end_date = dt.datetime(args.end_year, args.end_month, args.end_day, 23, 0)
height = args.height
lease_site = args.location

lease = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/lease_centroids.csv')
#lease = pd.read_csv('/Users/jameskim/Documents/rucool/Repositories/wind-science/files/lease_centroids.csv')
point = extract_coordinates(lease,lease_site)

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

# Filename and directory adjustment
safe_lease_site = lease_site.replace(" ", "")

directory_path = f'/www/web/rucool/windenergy/ru-wrf/wind_data/lease_centroids_csv/{safe_lease_site}/{height}'
#directory_path = f'/users/jameskim/Documents/rucool/Data/{safe_lease_site}/{height}'
csv_filename = f'{directory_path}/RUWRF_wind_data_{safe_lease_site}_{start_date.year}_{height}m.csv'

# Ensure directory exists
os.makedirs(directory_path, exist_ok=True)

# Save DataFrame to CSV
df.to_csv(csv_filename, index=False)
print(f'Data saved to {csv_filename}')