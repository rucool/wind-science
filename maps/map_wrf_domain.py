import numpy as np
import os
import glob
#from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#import geopandas as gpd
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.ticker as mticker
import functions.plotting as pf
import matplotlib.pyplot as plt
import cool_maps.plot as cplt

domainDir = '/home/wrfadmin/toolboxes/wind-science/files'
weaDir = '/home/coolgroup/bpu/mapdata/shapefiles/BOEM-Renewable-Energy-Shapefiles-current'
isoDir = '/home/coolgroup/bpu/mapdata/shapefiles/Independent_System_Operators_2019'
figDir = '/www/web/rucool/windenergy/ru-wrf/maps'
plot1 = True
plot3 = True
plot9 = True
plotWEA = True
plotISO = False

domain = [-75.5, -72, 38.4, 40.6]
coast='full'
if plot1:
    lon1=np.genfromtxt(os.path.join(domainDir,'XLONG_M_1km.csv'),delimiter=',')
    lat1=np.genfromtxt(os.path.join(domainDir,'XLAT_M_1km.csv'),delimiter=',')
    domain=[np.min(lon1)-.1,np.max(lon1)+.1,np.min(lat1)-.1,np.max(lat1)+.1]
if plot3:
    lon3=np.genfromtxt(os.path.join(domainDir,'XLONG_M_3km.csv'),delimiter=',')
    lat3=np.genfromtxt(os.path.join(domainDir,'XLAT_M_3km.csv'),delimiter=',')
    domain=[np.min(lon3)-.1,np.max(lon3)+.1,np.min(lat3)-.1,np.max(lat3)+.1]
    coast='high'
if plot9:
    lon9=np.genfromtxt(os.path.join(domainDir,'XLONG_M_9km.csv'),delimiter=',')
    lat9=np.genfromtxt(os.path.join(domainDir,'XLAT_M_9km.csv'),delimiter=',')
    domain=[np.min(lon9)-.1,np.max(lon9)+.1,np.min(lat9)-.1,np.max(lat9)+.1]
    coast='high'

if plotWEA:
    lease = glob.glob(os.path.join(weaDir,'BOEM_Wind_Leases_*.shp'))[0]
    plan = glob.glob(os.path.join(weaDir,'BOEM_Wind_Planning_Areas_*.shp'))[0]
    #leasing_areas = gpd.read_file(shape_file_lease)
    #leasing_areas = leasing_areas.to_crs(crs={'init': 'epsg:4326'})
    #planning_areas = gpd.read_file(shape_file_plan)
    #planning_areas = planning_areas.to_crs(crs={'init': 'epsg:4326'})

if plotISO:
    # shape_file_iso = '/Users/nazzaro/Documents/BPU/Independent_System_Operators/Independent_System_Operators.shp'
    # shape_file_iso = '/Users/nazzaro/Documents/BPU/iso_shape_test/iso_new.shp'
    # iso_areas = gpd.read_file(shape_file_iso)
    # iso_areas = iso_areas.to_crs(crs={'init': 'epsg:4326'})
    shape_file_iso = os.path.join(isoDir,'Independent_System_Operators.shp')
    iso_areas = iso_areas.to_crs(crs={'init': 'epsg:4326'})

#fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection=ccrs.PlateCarree()))
fig, ax = cplt.create(domain)
plt.rcParams.update({'font.size': 14})

# ax.set_extent(domain)  # [min lon, max lon, min lat, max lat]

#ax.set_xticks([-80,-76,-72,-68,-64,-60], crs=ccrs.PlateCarree())
#ax.set_yticks([32,35,38,41,44], crs=ccrs.PlateCarree())
#LAND = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='tan')

#state_lines = cfeature.NaturalEarthFeature(
#    category='cultural',
#    name='admin_1_states_provinces_lines',
#    scale='50m',
#    facecolor='none')


#ax.add_feature(LAND)#, edgecolor='black')
#ax.add_feature(LAND, edgecolor='black', linewidth=1.5, facecolor='none', zorder=7)
#ax.add_feature(cfeature.LAKES, facecolor='white')
#ax.add_feature(cfeature.BORDERS, color='black', linewidth=1.5, zorder=7)
#ax.add_feature(state_lines, zorder=7, edgecolor='black', linewidth=1.5)

figdesc=''
if plotWEA:
    #leasing_areas.plot(ax=ax, color='magenta', edgecolor='magenta')
    #planning_areas.plot(ax=ax, color='green', edgecolor='green')
    kwargs = dict()
    kwargs['edgecolor'] = 'magenta'
    kwargs['facecolor'] = 'magenta'
    pf.map_add_boem_outlines(ax, lease, **kwargs)
    kwargs = dict()
    kwargs['edgecolor'] = 'green'
    kwargs['facecolor'] = 'green'
    pf.map_add_boem_outlines(ax, plan, **kwargs)
    figdesc+='wea_'
if plotISO:
    # this uses geopandas - not supported in wind-science yet
    iso_areas[iso_areas['NAME']=='PJM INTERCONNECTION, LLC'].plot(ax=ax, color='#5D92B1', edgecolor='black', linewidth=.5, zorder=6)
    #iso_areas[iso_areas['NAME']=='MIDCONTINENT INDEPENDENT TRANSMISSION SYSTEM OPERATOR, INC..'].plot(ax=ax, color='#84BE6A', edgecolor='black', linewidth=.5, zorder=5)
    iso_areas[iso_areas['NAME']=='NEW YORK INDEPENDENT SYSTEM OPERATOR'].plot(ax=ax, color='#84BE6A', edgecolor='black', linewidth=.5, zorder=5)
    iso_areas[iso_areas['NAME']=='ISO NEW ENGLAND INC.'].plot(ax=ax, color='yellow', edgecolor='black', linewidth=.5, zorder=5)
    figdesc+='iso_'

lw=1.5
if plot3:
    ax.plot(lon3[0,:],lat3[0,:],color='red',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon3[-1,:],lat3[-1,:],color='red',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon3[:,0],lat3[:,0],color='red',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon3[:,-1],lat3[:,-1],color='red',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
if plot9:
    ax.plot(lon9[0,:],lat9[0,:],color='blue',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon9[-1,:],lat9[-1,:],color='blue',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon9[:,0],lat9[:,0],color='blue',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon9[:,-1],lat9[:,-1],color='blue',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
if plot1:
    ax.plot(lon1[0,:],lat1[0,:],color='navy',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon1[-1,:],lat1[-1,:],color='navy',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon1[:,0],lat1[:,0],color='navy',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
    ax.plot(lon1[:,-1],lat1[:,-1],color='navy',linewidth=lw, transform=ccrs.PlateCarree(),zorder=20)
if plot9:
    ax.text(lon9[0,-1]-.1,lat9[0,-1]+.3,'RUWRF 9km',color='blue',horizontalalignment='right',verticalalignment='bottom',fontsize=14, transform=ccrs.PlateCarree(),zorder=20)
    figdesc+='9km_'
if plot3:
    ax.text(lon3[0,-1]-.1,lat3[0,-1]+.1,'RUWRF 3km',color='red',horizontalalignment='right',verticalalignment='bottom',fontsize=14, transform=ccrs.PlateCarree(),zorder=20)
    figdesc+='3km_'
if plot1:
    ax.text(lon1[0,-1]-.1,lat1[0,-1]+.1,'RUWRF 1km',color='navy',horizontalalignment='right',verticalalignment='bottom',fontsize=6, transform=ccrs.PlateCarree(),zorder=20)
    figdesc+='1km_'

plt.savefig(os.path.join(figDir,'RUWRF_domain_'+figdesc+pd.to_datetime('now').strftime('%Y%m%d')+'.png'),dpi=300)
