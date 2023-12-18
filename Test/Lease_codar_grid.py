import os
import pyproj
import xarray as xr
import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from shapely.geometry import Point
import functions.common as cf
import functions.plotting as pf
import cool_maps.plot as cplt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import math

plt.rcParams.update({'font.size': 12})

# Function to load latitude and longitude points from a file
def load_points(file_path):
    data = np.loadtxt(file_path, delimiter=' ')
    lon_p = data[:, 1]
    lat_p = data[:, 0]
    return lat_p, lon_p

# Function to generate range rings
def generate_range_rings(rc_size, num_of_range_cells, angular_coverage, start_angle, radial_spacing, lat, lon):
    angular_bins = int(angular_coverage / radial_spacing)
    angular_coverage = angular_coverage * math.pi / 180

    r_ring_x = np.zeros((angular_bins + 1, num_of_range_cells))
    r_ring_y = np.zeros((angular_bins + 1, num_of_range_cells))

    for i in range(angular_bins + 1):
        for r in range(num_of_range_cells):
            clat = lat + ((r * rc_size - rc_size/2) * math.sin(((start_angle + 82.5) * math.pi/180) + i * angular_coverage / angular_bins)) / 111.12
            clon = lon - ((r * rc_size - rc_size/2) * math.cos(((start_angle + 82.5) * math.pi/180) + i * angular_coverage / angular_bins) / (111.12 * math.cos(lat * math.pi/180)))
            r_ring_x[i, r] = clon
            r_ring_y[i, r] = clat

    r_line_x = np.zeros((num_of_range_cells, angular_bins + 1))
    r_line_y = np.zeros((num_of_range_cells, angular_bins + 1))

    for r in range(num_of_range_cells):
        for i in range(angular_bins + 1):
            clat2 = lat + ((r * rc_size - rc_size/2) * math.sin(((start_angle + 82.5) * math.pi/180) + i * angular_coverage / angular_bins)) / 111.12
            clon2 = lon - ((r * rc_size - rc_size/2) * math.cos(((start_angle + 82.5) * math.pi/180) + i * angular_coverage / angular_bins) / (111.12 * math.cos(lat * math.pi/180)))
            r_line_x[r, i] = clon2
            r_line_y[r, i] = clat2

    return r_ring_x, r_ring_y, r_line_x, r_line_y

# Function to plot radar sites
def plot_radar_site(ax, site_name, lat, lon, rc_size, num_of_range_cells, angular_coverage, start_angle, radial_spacing):
    r_ring_x, r_ring_y, r_line_x, r_line_y = generate_range_rings(rc_size, num_of_range_cells, angular_coverage, start_angle, radial_spacing, lat, lon)

    for i in range(r_ring_x.shape[0]):
        ax.plot(r_ring_x[i], r_ring_y[i], color='black', linewidth=0.8, alpha=0.50, transform=ccrs.PlateCarree())

    for r in range(r_line_x.shape[0]):
        ax.plot(r_line_x[r], r_line_y[r], color='red', linewidth=0.8, alpha=0.50, transform=ccrs.PlateCarree())

    ax.scatter(lon, lat, color='yellow', marker='^', s=100, zorder=1000, transform=ccrs.PlateCarree(), label=site_name)
    ax.text(lon - 0.02, lat + 0.02, site_name, transform=ccrs.PlateCarree(), fontsize=14, fontweight='bold', zorder=2000,
            verticalalignment='bottom', horizontalalignment='right')

# Function to create legend
def create_legend(ax):
    legend = ax.legend(loc='upper right', frameon=True, edgecolor='black', fontsize='medium', markerscale=2,)
    for handle in legend.legendHandles:
        handle.set_sizes([40])
    return legend

# Function to create the map
def create_map(extent):
    fig, ax = cplt.create(extent, proj=ccrs.Mercator(), bathymetry=False, isobaths=None, figsize=(16, 9))
    return fig, ax

# Main code
# Martha's Vineyard Site
extent_marthas = [-71.1, -69.95, 40.7, 41.5]  # Evangrid
fig_marthas, ax_marthas = create_map(extent_marthas)

lease_marthas = '/Users/jameskim/documents/rucool/WRF/Wind_Lease_Outlines_2_2023.shp'
bwargs_marthas = dict(edgecolor='white', zorder=10, linewidth=1.5)
geoms_marthas = pf.map_add_boem_outlines(ax_marthas, lease_marthas, **bwargs_marthas)

df_cod = pd.read_csv('/users/jameskim/Documents/rucool/CODAR/data/maracoos_codar_sites_2023.csv')
lat_marthas_codar = df_cod.loc[2, ' Latitude']
lon_marthas_codar = df_cod.loc[2, ' Longitude']
label_marthas_codar = 'MVCO'

ax_marthas.scatter(lon_marthas_codar, lat_marthas_codar, color='yellow', marker='^', s=100, zorder=1000,
                   transform=ccrs.PlateCarree())
ax_marthas.text(lon_marthas_codar - 0.02, lat_marthas_codar + 0.02, label_marthas_codar,
                transform=ccrs.PlateCarree(), fontsize=14, fontweight='bold', zorder=2000,
                verticalalignment='bottom', horizontalalignment='right')

rc_size_marthas = 6
num_of_range_cells_marthas = 30
angular_coverage_marthas = 140
start_angle_marthas = 115
radial_spacing_marthas = 5

plot_radar_site(ax_marthas, 'Martha\'s Vineyard', lat_marthas_codar, lon_marthas_codar, rc_size_marthas,
                 num_of_range_cells_marthas, angular_coverage_marthas, start_angle_marthas, radial_spacing_marthas)

create_legend(ax_marthas)

ax_marthas.add_feature(cfeature.COASTLINE, edgecolor='black')
ax_marthas.add_feature(cfeature.BORDERS, linestyle=':')
ax_marthas.set_title('Martha\'s Vineyard Site', fontsize=20, fontweight='bold')

plt.show()

############################################################################################################
############## Nantucket Site ##############################################################################
############################################################################################################
extent_nantucket = [-72.0,-69.5,40.3,41.5] 
fig_nantucket, ax_nantucket = create_map(extent_nantucket)

lease_nantucket = '/Users/jameskim/documents/rucool/WRF/Wind_Lease_Outlines_2_2023.shp'
bwargs_nantucket = dict(edgecolor='white', zorder=10, linewidth=1.5)
geoms_nantucket = pf.map_add_boem_outlines(ax_nantucket, lease_nantucket, **bwargs_nantucket)

df_cod = pd.read_csv('/users/jameskim/Documents/rucool/CODAR/data/maracoos_codar_sites_2023.csv')
lat_nantucket_codar = df_cod.loc[1, ' Latitude']
lon_nantucket_codar = df_cod.loc[1, ' Longitude']
label_nantucket_codar = 'Nantucket'

ax_nantucket.scatter(lon_nantucket_codar, lat_nantucket_codar, color='yellow', marker='^', s=100, zorder=1000,
                     transform=ccrs.PlateCarree())
ax_nantucket.text(lon_nantucket_codar - 0.02, lat_nantucket_codar + 0.02, label_nantucket_codar,
                  transform=ccrs.PlateCarree(), fontsize=14, fontweight='bold', zorder=2000,
                  verticalalignment='bottom', horizontalalignment='right')

rc_size_nantucket = 6
num_of_range_cells_nantucket = 30
angular_coverage_nantucket = 140
start_angle_nantucket = 115
radial_spacing_nantucket = 5

plot_radar_site(ax_nantucket, 'Nantucket', lat_nantucket_codar, lon_nantucket_codar, rc_size_nantucket,
                 num_of_range_cells_nantucket, angular_coverage_nantucket, start_angle_nantucket,
                 radial_spacing_nantucket)

create_legend(ax_nantucket)

ax_nantucket.add_feature(cfeature.COASTLINE, edgecolor='black')
ax_nantucket.add_feature(cfeature.BORDERS, linestyle=':')
ax_nantucket.set_title('Nantucket Site', fontsize=20, fontweight='bold')

plt.show()