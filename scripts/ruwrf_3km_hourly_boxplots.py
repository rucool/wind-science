from datetime import timedelta, datetime as dt
import os
import xarray as xr
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import functions.common as cf
import calendar
import argparse
import sys
import pytz
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified

def get_month_dates(month=None, season=None, start_year=None, end_year=None):
    if month and season:
        raise ValueError("Please provide either 'month' or 'season', not both.")

    results = []

    if month:
        # Convert month name to its numerical representation
        month_number = dt.strptime(month, '%B').month

        for year in range(start_year, end_year + 1):
            # Calculate the first day of the month
            start_date = dt(year, month_number, 1)

            # Calculate the last day of the month
            if month_number == 12:
                end_date = dt(year + 1, 1, 1) - timedelta(days=1)
            else:
                end_date = dt(year, month_number + 1, 1) - timedelta(days=1)

            # Format dates as strings
            start_date_str = start_date.strftime('%Y-%m-%d')
            end_date_str = end_date.strftime('%Y-%m-%d')

            # Create a dictionary with column names and values
            result_dict = {
                'start_date': start_date_str,
                'end_date': end_date_str,
                'month_year': f"{month} {year}"
            }

            results.append(result_dict)

    elif season:
        # Define seasons and their corresponding months
        seasons = {
            'spring': [3, 4, 5],
            'summer': [6, 7, 8],
            'fall': [9, 10, 11],
            'winter': [12, 1, 2]
        }

        for year in range(start_year, end_year + 1):
            # Get the months corresponding to the season
            season_months = seasons[season.lower()]

            # Calculate the first day of the season
            start_date = dt(year - 1, 12, 1) if 'winter' in season.lower() else dt(year, min(season_months), 1)

            # Calculate the last day of the season
            if max(season_months) == 12:
                end_date = dt(year, 2, 28)  # Assuming February has 28 days
            else:
                end_date = dt(year, max(season_months) + 1, 1) - timedelta(days=1)

            # Format dates as strings
            start_date_str = start_date.strftime('%Y-%m-%d')
            end_date_str = end_date.strftime('%Y-%m-%d')

            # Create a dictionary with column names and values
            result_dict = {
                'start_date': start_date_str,
                'end_date': end_date_str,
                'season_year': f"{season.capitalize()} {year}"
            }

            results.append(result_dict)

    return results

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

def process_wind_data(ds, point, start_year, end_year, season, month, power_curve):
    # Initialize the overall DataFrame
    ws_df_all = pd.DataFrame()

    # Iterate over each year in the specified range
    for year in range(start_year, end_year + 1):
        # Get month dates based on the specified season and month
        dates = get_month_dates(season=season, month=month, start_year=year, end_year=year)

        # Create an empty DataFrame for the current year
        ws_df_year = pd.DataFrame()

        # Iterate over each date range in the current year
        for date_range in dates:
            start_date = dt.strptime(date_range['start_date'], '%Y-%m-%d')
            end_date = dt.strptime(date_range['end_date'], '%Y-%m-%d')

            # Slice the dataset for the current date range
            ds_subset = ds.sel(time=slice(start_date, end_date))

            # Find the closest grid point to the specified location
            a = abs(ds_subset['XLAT'] - point[0]) + abs(ds_subset['XLONG'] - point[1])
            i, j = np.unravel_index(a.argmin(), a.shape)

            # Extract U and V components at the specified height
            u = ds_subset.sel(height=160)['U'][:, i, j]
            v = ds_subset.sel(height=160)['V'][:, i, j]

            # Calculate wind speed
            ws = cf.wind_uv_to_spd(u, v)

            # Create a DataFrame for the current date range
            ws_df_range = ws.to_dataframe('windspeed')
            ws_df_range['Power'] = np.interp(ws_df_range['windspeed'], power_curve['Wind Speed'], power_curve['Power'])

            # Concatenate the results to the DataFrame for the current year
            ws_df_year = pd.concat([ws_df_year, ws_df_range])

        # Concatenate the results for the current year to the overall DataFrame
        ws_df_all = pd.concat([ws_df_all, ws_df_year])

    # Extract hour from the index
    ws_df_all['hour'] = ws_df_all.index.hour

    # Drop the first row
    ws_df_all = ws_df_all.drop(ws_df_all.index[0])

    # Reset the index to avoid issues when running the code multiple times
    ws_df_all = ws_df_all.reset_index(drop=True)

    # Change all occurrences of hour 0 to 24
    ws_df_all['hour'] = ws_df_all['hour'].replace(0, 24)

    return ws_df_all

def plot_and_save_power_variation(ws_df_all, lease_site, start_year, end_year, season, month, output_directory):
    
    medianprops = dict(linestyle="-", color='black', linewidth=2)
    meanpointprops = dict(marker='D', markeredgecolor='black', markerfacecolor='black')

    plt.figure(figsize=(10, 6))
    box = plt.boxplot([group['Power'] for name, group in ws_df_all.groupby('hour')],
                     showfliers=True, patch_artist=False, showmeans=True,
                     medianprops=medianprops, meanprops=meanpointprops)

    plt.xlabel('Hour of the Day (EST)', fontsize=14)
    plt.ylabel('Power (kw)', fontsize=14)

    if start_year == end_year and season is None:
        plt.title(f'Hourly Power (160m) Variation at {lease_site} ({month} {start_year})', fontsize=16, fontweight='bold')
    elif start_year != end_year and season is None:
        plt.title(f'Hourly Power (160m) Variation at {lease_site} ({month} {start_year}-{end_year})', fontsize=16, fontweight='bold')

    if start_year == end_year and month is None:
        plt.title(f'Hourly Power (160m) Variation at {lease_site} ({season} {start_year})', fontsize=16, fontweight='bold')
    elif start_year != end_year and month is None:
        plt.title(f'Hourly Power (160m) Variation at {lease_site} ({season} {start_year}-{end_year})', fontsize=16, fontweight='bold')

    # Define a fixed offset of -4 hours (Eastern Daylight Time)
    fixed_offset = pytz.FixedOffset(-240)  # -240 minutes = -4 hours

    # Assuming you have the labels generated before the timezone conversion
    labels = [str(((i - 1) % 12) + 1) + (' AM' if i < 12 else ' PM') for i in range(1, 25, 3)]

    # Convert labels to datetime objects in GMT
    gmt_labels = [dt.strptime(label, '%I %p').replace(tzinfo=pytz.timezone('GMT')) for label in labels]

    # Convert GMT labels to Eastern Time with a fixed offset
    est_labels = [dt.astimezone(fixed_offset).strftime('%-I %p') for dt in gmt_labels]

    # Convert to coordinate time with fixed offset
    plt.xticks(ticks=np.arange(1, 25, 3), labels=est_labels)

    # Customizing gridlines
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # Adding a background grid to the plot
    plt.gca().set_facecolor('#f7f7f7')

    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()

    # Set y-axis limits with a little padding
    plt.ylim(-200, 15500)

    # Set y-axis ticks manually
    plt.yticks([0, 3000, 6000, 9000, 12000, 15000])

    # Legend
    median_legend = mlines.Line2D([], [], color='black', label='Median')
    mean_legend = mlines.Line2D([], [], marker='D', markerfacecolor='black', markeredgecolor='black', markersize=10,
                               label='Mean', linestyle='None')

    plt.legend(handles=[median_legend, mean_legend], loc='upper right')

    # Save the plot to the specified directory
    # Create the output directory based on arguments
    if season:
        output_directory = os.path.join(output_directory, f'power/{lease_site}/season')
    elif start_year == end_year and month:
        output_directory = os.path.join(output_directory, f'power/{lease_site}/months')
    elif start_year != end_year and month:
        output_directory = os.path.join(output_directory, f'power/{lease_site}/monthly')

 # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Create the filename based on the arguments
    if start_year == end_year:
        filename = f'3km_power_hourlybox_{lease_site}_{start_year}'
        filename += f'_{month}.png' if month else f'_{season}.png'
    else:
        filename = f'3km_power_hourlybox_{lease_site}_{start_year}_{end_year}'
        filename += f'_{month}.png' if month else f'_{season}.png'

    filepath = os.path.join(output_directory, filename)

    # Save the plot to the specified directory
    plt.savefig(filepath)
 
def plot_and_save_wind_speed_variation(ws_df_all, lease_site, start_year, end_year, season, month, output_directory):
    
    plt.clf()

    medianprops = dict(linestyle="-", color='black', linewidth=2)
    meanpointprops = dict(marker='D', markeredgecolor='black', markerfacecolor='black')

    plt.figure(figsize=(10, 6))
    box = plt.boxplot([group['windspeed'] for name, group in ws_df_all.groupby('hour')],
                      showfliers=True, patch_artist=False, showmeans=True,
                      medianprops=medianprops, meanprops=meanpointprops)

    plt.xlabel('Hour of the Day (EST)', fontsize=14)
    plt.ylabel('Wind Speed (m/s)', fontsize=14)

    if start_year == end_year and season is None:
        plt.title(f'Hourly Wind Speed (160m) Variation at {lease_site} ({month} {start_year})', fontsize=16,
                  fontweight='bold')
    elif start_year != end_year and season is None:
        plt.title(f'Hourly Wind Speed (160m) Variation at {lease_site} ({month} {start_year}-{end_year})', fontsize=16,
                  fontweight='bold')

    if start_year == end_year and month is None:
        plt.title(f'Hourly Wind Speed (160m) Variation at {lease_site} ({season} {start_year})', fontsize=16,
                  fontweight='bold')
    elif start_year != end_year and month is None:
        plt.title(f'Hourly Wind Speed (160m) Variation at {lease_site} ({season} {start_year}-{end_year})', fontsize=16,
                  fontweight='bold')

    # Define a fixed offset of -4 hours (Eastern Daylight Time)
    fixed_offset = pytz.FixedOffset(-240)  # -240 minutes = -4 hours

    # Assuming you have the labels generated before the timezone conversion
    labels = [str(((i - 1) % 12) + 1) + (' AM' if i < 12 else ' PM') for i in range(1, 25, 3)]

    # Convert labels to datetime objects in GMT
    gmt_labels = [dt.strptime(label, '%I %p').replace(tzinfo=pytz.timezone('GMT')) for label in labels]

    # Convert GMT labels to Eastern Time with a fixed offset
    est_labels = [dt.astimezone(fixed_offset).strftime('%-I %p') for dt in gmt_labels]

    # Convert to coordinate time with fixed offset
    plt.xticks(ticks=np.arange(1, 25, 3), labels=est_labels)

    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    plt.gca().set_facecolor('#f7f7f7')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    plt.ylim(0, 30)

    median_legend = mlines.Line2D([], [], color='black', label='Median')
    mean_legend = mlines.Line2D([], [], marker='D', markerfacecolor='black', markeredgecolor='black', markersize=10,
                                label='Mean', linestyle='None')

    plt.legend(handles=[median_legend, mean_legend], loc='upper right')

    # Save the plot to the specified directory
    # Create the output directory based on arguments
    if season:
        output_directory = os.path.join(output_directory, f'windspeed/{lease_site}/season')
    elif start_year == end_year and month:
        output_directory = os.path.join(output_directory, f'windspeed/{lease_site}/months')
    elif start_year != end_year and month:
        output_directory = os.path.join(output_directory, f'windspeed/{lease_site}/monthly')

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)


    # Create the filename based on the arguments
    if start_year == end_year:
        filename = f'3km_windspeed_hourlybox_{lease_site}_{start_year}'
        filename += f'_{month}.png' if month else f'_{season}.png'
    else:
        filename = f'3km_windspeed_hourlybox_{lease_site}_{start_year}_{end_year}'
        filename += f'_{month}.png' if month else f'_{season}.png'

    filepath = os.path.join(output_directory, filename)

    # Save the plot to the specified directory
    plt.savefig(filepath)




################################ Main Code ################################

def main(args):
    lease = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/lease_centroids.csv') #
    power_curve = pd.read_csv('/home/wrfadmin/toolboxes/wind-science/files/wrf_lw15mw_power.csv') 
    # /Users/jameskim/Documents/rucool/Repositories/wind-science/files/lease_centroids.csv
    # /Users/jameskim/Documents/rucool/Repositories/wind-science/files/wrf_lw15mw_power.csv



    lease_site = args.location
    season = args.season
    month = args.month
    start_year = args.start_year
    end_year = args.end_year

    # import WRF data link
    mlink = 'http://tds.marine.rutgers.edu/thredds/dodsC/cool/ruwrf/wrf_4_1_3km_processed/WRF_4.1_3km_Processed_Dataset_Best'
    ds = xr.open_dataset(mlink)

    # Create empty DataFrames to store windspeed and power data
    ws_df_all = pd.DataFrame()
    power_df_all = pd.DataFrame()

    # HERE WE GO WOAHHHH
    
    point = extract_coordinates(lease,lease_site)
    df = process_wind_data(ds, point, start_year, end_year, season, month, power_curve)
    plot_and_save_power_variation(df, lease_site, start_year, end_year, season, month, output_directory = '/www/web/rucool/windenergy/ru-wrf/images/hourly_boxplots') 
    plot_and_save_wind_speed_variation(df, lease_site, start_year, end_year, season, month, output_directory = '/www/web/rucool/windenergy/ru-wrf/images/hourly_boxplots') 
    
    
    #  /www/web/rucool/windenergy/ru-wrf/images/hourly_boxplots
    #   /users/jameskim/documents/test_figures
if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=main.__doc__,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-s', '--start',
                            dest='start_year',
                            default='2022',
                            type=int,
                            help='Date format in YYYY')

    arg_parser.add_argument('-e', '--end',
                            dest='end_year',
                            default='2022',
                            type=int,
                            help='Date in format YYYY')

    arg_parser.add_argument('-location',
                            dest='location',
                            default='OCS-A 0499',
                            type=str,
                            help='Lease area point at which to grab wind speeds to generate wind rose. Valid OCS lease '
                                 'codes can be found in ./files/lease_centroids.csv. Example: "OCS-A0499".')
    
    arg_parser.add_argument('-sn','--season',
                            dest='season',
                            default = None,
                            type=str,
                            help='Enter name of Season')
    
    arg_parser.add_argument('-m','--month',
                            dest='month',
                            default = None,
                            type=str,
                            help='Enter name of Season')
    
    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))