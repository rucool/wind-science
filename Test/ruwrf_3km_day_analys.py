# Author: James Kim
# Date created 11/10/23 | Last Modified:  

# Generate baseline products for hourly wind and power on a daily timescale from 3km ruwrf

import datetime as dt
import os
import xarray as xr
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import functions.common as cf
import calendar
import pytz
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified

#arguments

#point = [39.17745, -74.18033]
lease_site = "OCS-A0000"
lat = 39.4
lon = -73.8
start_date = dt.datetime(2020, 8, 1, 0, 0)  
end_date = dt.datetime(2020, 8, 31, 23, 0)

# grab the location of the specified site
# need to update this #
point = [lat,lon]



# import WRF data link 
mlink = 'http://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
ds = xr.open_dataset(mlink)
ds = ds.sel(time=slice(start_date, end_date))


# Find the nearest latitude and longitude in the dataset

 # calculate the sum of the absolute value distance between the model location and buoy location
a = abs(ds['XLAT'] - point[0]) + abs(ds['XLONG'] - point[1])

# find the indices of the minimum value in the array calculated above
i, j = np.unravel_index(a.argmin(), a.shape)

# get u and v component, calculate windspeed and estimated wind power
u = ds.sel(height=160)['U'][:, i, j]
v = ds.sel(height=160)['V'][:, i, j]
ws = cf.wind_uv_to_spd(u,v)

ws_df = ws.to_dataframe('windspeed')

# create Power
power_curve =pd.read_csv('/Users/jameskim/Documents/rucool/Repositories/wind-science/files/wrf_lw15mw_power.csv')

# Interpolate power values and add a new column 'Power' to ws_df
ws_df['Power'] = np.interp(ws_df['windspeed'], power_curve['Wind Speed'], power_curve['Power'])

# Extract hour from the index of ws_df and convert to EST
est = pytz.timezone('US/Eastern')  # Define the Eastern Timezone
ws_df.index = ws_df.index.tz_localize('GMT').tz_convert(est)  # Convert time to EST

# Extract hour from the index of ws_df
ws_df['hour'] = ws_df.index.hour


# customize the boxplot elements


# Clear the current figure

plt.clf()

medianprops = dict(linestyle = "-", color='black',linewidth  = 2)
meanpointprops = dict(marker='D', markeredgecolor='black', markerfacecolor='black')

plt.figure(figsize=(10, 6))
box = plt.boxplot([group['windspeed'] for name, group in ws_df.groupby('hour')], showfliers=True, patch_artist=False,showmeans=True,medianprops=medianprops,meanprops = meanpointprops)

plt.xlabel('Hour of the Day (EST)', fontsize=12)
plt.ylabel('Wind Speed (m/s)', fontsize=12)
plt.title(f'Hourly Wind Speed Variation at OCS-0000 (Summer 2020)', fontsize=14)

# convert to coordinal time
plt.xticks(ticks=np.arange(1, 25, 3), labels=[str((i % 12) if i % 12 else 12) + (' PM' if i >= 12 else ' AM') for i in range(0, 24, 3)])

# Customizing gridlines
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

# Adding a background grid to the plot
plt.gca().set_facecolor('#f7f7f7')

plt.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()

# Setting y-axis limits
plt.ylim(0, 30)

##### legend ####

# Create legends for the median and mean without the additional line
median_legend = mlines.Line2D([], [], color='black', label='Median')
mean_legend = mlines.Line2D([], [], marker='D', markerfacecolor='black', markeredgecolor='black', markersize=10, label='Mean', linestyle='None')

plt.legend(handles=[median_legend, mean_legend], loc='upper right')

plt.show()




# Clear the current figure
plt.clf()

medianprops = dict(linestyle="-", color='black', linewidth=2)
meanpointprops = dict(marker='D', markeredgecolor='black', markerfacecolor='black')

plt.figure(figsize=(10, 6))
box = plt.boxplot([group['Power'] for name, group in ws_df.groupby('hour')],
                 showfliers=True, patch_artist=False, showmeans=True,
                 medianprops=medianprops, meanprops=meanpointprops)

plt.xlabel('Hour of the Day (EST)', fontsize=14)
plt.ylabel('Power (kw)', fontsize=14)
plt.title(f'Hourly Power Variation at OCS-0000 (Summer 2020)', fontsize=14, fontweight = "bold")

# Convert to coordinate time
plt.xticks(ticks=np.arange(1, 25, 3),
           labels=[str((i % 12) if i % 12 else 12) + (' PM' if i >= 12 else ' AM') for i in range(0, 24, 3)])

# Customizing gridlines
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

# Adding a background grid to the plot
plt.gca().set_facecolor('#f7f7f7')

plt.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()

# Set y-axis limits with a little padding
plt.ylim(-200, 15500)

# Set y-axis ticks manually
plt.yticks([0,3000, 6000, 9000, 12000, 15000])

# Legend
median_legend = mlines.Line2D([], [], color='black', label='Median')
mean_legend = mlines.Line2D([], [], marker='D', markerfacecolor='black', markeredgecolor='black', markersize=10, label='Mean', linestyle='None')

plt.legend(handles=[median_legend, mean_legend], loc='upper right')

plt.show()















