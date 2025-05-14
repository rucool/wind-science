#!/usr/bin/env python

"""
Author: Lori Garzio on 9/26/2024
Last modified: 5/14/2025
Plot surface maps of Summer 2022 averaged datasets for variables extracted from two versions of RU-WRF model output:
a version with simulated turbines compared to a control version with no turbines. Also plot the calculated differences
(turbines minus control).
"""

import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import cmocean as cmo
from math import ceil
import cool_maps.plot as cplt
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def add_contours(ax, londata, latdata, vardata, clist, label_format=None, color='black', linewidth=.5, linestyle='-'):
    """
    Adds black contour lines with labels to a cartopy map object
    :param ax: plotting axis object
    :param londata: longitude data
    :param latdata: latitude data
    :param vardata: variable data
    :param clist: list of contour levels
    :param label_format: optional format for contour labels (e.g. '%.1f')
    """
    label_format = label_format or '%d'

    CS = ax.contour(londata, latdata, vardata, clist, colors=color, linewidths=linewidth,
                    linestyle=linestyle, transform=ccrs.PlateCarree())
    ax.clabel(CS, inline=True, fontsize=10.5, fmt=label_format)


def subset_grid(ext, data, lon_name, lat_name):
    """
    This function subsets the model output to the defined boundaries
    ext: list of boundaries [lon_min, lon_max, lat_min, lat_max]
    data: xarray dataset array
    lon_name: name of the longitude coordinate
    lat_name: name of the latitude coordinate
    """
    lonx = data[lon_name]
    laty = data[lat_name]

    lon_idx = np.logical_and(lonx > ext[0], lonx < ext[1])
    lat_idx = np.logical_and(laty > ext[2], laty < ext[3])

    # find i and j indices of lon/lat in boundaries
    ind = np.where(np.logical_and(lat_idx, lon_idx))

    # subset data from min i,j corner to max i,j corner
    # there will be some points outside of defined boundaries because grid is not rectangular
    data_sub = np.squeeze(data)[range(np.min(ind[0]), np.max(ind[0]) + 1), range(np.min(ind[1]), np.max(ind[1]) + 1)]
    lon = data_sub[lon_name]
    lat = data_sub[lat_name]

    return data_sub, lon, lat


def map_add_boem_outlines(ax, shpfile, edgecolor=None, facecolor=None, projection=None, zorder=None, alpha=None,
                          linewidth=1):
    edgecolor = edgecolor or 'black'
    facecolor = facecolor or 'none'
    projection = projection or ccrs.PlateCarree()
    zorder = zorder or 1
    alpha = alpha or 1

    shape_feature = cfeature.ShapelyFeature(
        Reader(shpfile).geometries(),
        projection,
        edgecolor=edgecolor,
        facecolor=facecolor,
        linewidth=linewidth
        )

    ax.add_feature(shape_feature, zorder=zorder, alpha=alpha)

    return shape_feature._geoms


def plot_pcolormesh(fig, ax, lon_data, lat_data, var_data, cmap='jet', clab=None, var_lims=None,
                    norm_clevs=None, extend='neither', cbar_ticks=None, shading='nearest'):
    """
    Create a pseudocolor plot
    :param fig: figure object
    :param ax: plotting axis object
    :param lon_data: longitude data
    :param lat_data: latitude data
    :param var_data: variable data
    :param cmap: optional color map, default is 'jet'
    :param clab: optional colorbar label
    :param var_lims: optional, [minimum, maximum] values for plotting (for fixed colorbar)
    :param norm_clevs: optional normalized levels
    :param extend: optional extend the colorbar, default is 'neither'
    """

    plt.subplots_adjust(right=0.88)
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    fig.add_axes(cax)

    if var_lims:
        h = ax.pcolormesh(lon_data, lat_data, var_data, vmin=var_lims[0], vmax=var_lims[1],  cmap=cmap,
                          transform=ccrs.PlateCarree())
    elif norm_clevs:
        h = ax.pcolormesh(lon_data, lat_data, var_data, cmap=cmap, norm=norm_clevs,
                          transform=ccrs.PlateCarree(), shading=shading)
    else:
        h = ax.pcolormesh(lon_data, lat_data, var_data, cmap=cmap, transform=ccrs.PlateCarree())

    cb = plt.colorbar(h, cax=cax, extend=extend)
    if clab:
        cb.set_label(label=clab, fontsize=14)
    if cbar_ticks:
        cb.set_ticks(cbar_ticks)


def main(variable, h, coastline, save_dir):
    extent = [-75.1, -72.2, 38.1, 40.5]

    f = f'/Users/garzio/Documents/rucool/Miles/RMI/2024_Ocean_Mixing/data/1km_wf2km_nyb/ruwrf_1km_wf2km_nyb_{variable}_avg_20220601_20220831.nc'
    fc = f'/Users/garzio/Documents/rucool/Miles/RMI/2024_Ocean_Mixing/data/1km_ctrl/ruwrf_1km_ctrl_{variable}_avg_20220601_20220831.nc'

    ds = xr.open_dataset(f)
    dsc = xr.open_dataset(fc)

    try:
        ds = ds.sel(height=h)
        dsc = dsc.sel(height=h)
    except KeyError:
        print('')

    varnameavg = f'{variable}_avg'

    configs = dict(
        ws=dict(cmap=plt.get_cmap('BuPu'), diffcmap=plt.get_cmap('RdBu_r'), diffmask=[-0.1, 0.1], difflims=[-2, 2], diffbins=40),
        ws10=dict(cmap=plt.get_cmap('BuPu'), diffcmap=plt.get_cmap('RdBu_r'), diffmask=[-0.1, 0.1], difflims=[-2, 2], diffbins=40),
        #UST=dict(cmap=plt.get_cmap('YlGnBu'), diffcmap=plt.get_cmap('Spectral_r'), diffmask=[-0.01, 0.01], lims=[0.1, 0.3], bins=20, difflims=[-0.1, 0.1], diffbins=60),
        UST=dict(cmap=plt.get_cmap('YlGnBu'), diffcmap=plt.get_cmap('RdBu_r'), diffmask=[-0.001, 0.001], lims=[0.1, 0.3], bins=20, difflims=[-0.04, 0.04], diffbins=60),
        Q2=dict(cmap=cmo.cm.rain, diffcmap=plt.get_cmap('BrBG'), diffmask=[-0.000025, 0.000025], lims=[0, 0.02], bins=20, difflims=[-0.00015, 0.00015], diffbins=20),
        T2=dict(cmap=cmo.cm.thermal, diffmask=[-0.01, 0.01], difflims=[-0.2, 0.2], diffbins=40),
        TEMP=dict(cmap=cmo.cm.thermal, diffmask=[-0.01, 0.01], difflims=[-0.2, 0.2], diffbins=40),
        TKE_PBL=dict(cmap=cmo.cm.tempo, diffmask=[-0.01, 0.01], lims=[0, 0.2], bins=20, difflims=[-0.12, 0.12], diffbins=40),
        HFX=dict(cmap=cmo.cm.thermal, diffmask=[-0.1, 0.1], lims=[-10, 10], bins=100, difflims=[-1.5, 1.5], diffbins=40),
        TSK=dict(cmap=cmo.cm.thermal, diffmask=[-0.01, 0.01], difflims=[-0.2, 0.2], diffbins=40)
    )

    vconfigs = configs[v]
    vsub, lon, lat = subset_grid(extent, ds[varnameavg], 'XLONG', 'XLAT')
    vctrlsub, lon, lat = subset_grid(extent, dsc[varnameavg], 'XLONG', 'XLAT')

    #if v in ['TEMP', 'T2', 'TSK']:
    #    ln = vsub.attrs['long_name']
    #    vsub = vsub - 273.15  # convert from K to degrees C
    #    vsub.attrs['units'] = 'deg_C'
    #    vsub.attrs['long_name'] = ln
    #    vctrlsub = vctrlsub - 273.15  # convert from K to degrees C
    #    vctrlsub.attrs['units'] = 'deg_C'
    #    vctrlsub.attrs['long_name'] = ln

    diff = vsub - vctrlsub
    print(f'diff range: [{str(np.round(np.nanmin(diff), 3))}, {str(np.round(np.nanmax(diff), 3))}]')
    diff.attrs['units'] = vsub.attrs['units']
    diff.attrs['long_name'] = f'{vsub.attrs["long_name"]} Difference'
    try:
        masked_diff = np.ma.masked_inside(diff, vconfigs['diffmask'][0], vconfigs['diffmask'][1])
    except KeyError:
        masked_diff = diff

    # define color map and label
    try:
        cmin = vconfigs['lims'][0]
        cmax = vconfigs['lims'][1]
        bins = vconfigs['bins']
    except KeyError:
        cmin = np.floor(np.nanmin([vsub, vctrlsub]))
        cmax = np.ceil(np.nanmax([vsub, vctrlsub]))
        bins = (cmax - cmin) * 5

    cmap = vconfigs['cmap']
    levels = MaxNLocator(nbins=bins).tick_values(cmin, cmax)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    if vsub.units == 'm s-1':
        unitslab = '($\\mathrm{{m}}\\ \\mathrm{{s}}^{{-1}}$)'
    elif vsub.units == 'kg kg-1':
        unitslab = '($\\mathrm{{kg}}\\ \\mathrm{{kg}}^{{-1}}$)'
    else:
        unitslab = f'({vsub.units})'
    color_label = f'{vsub.long_name} {unitslab}'
######################################################################################################################
    # plot average: wind farm
    kwargs = dict()
    kwargs['figsize'] = (8, 8)
    kwargs['oceancolor'] = 'none'
    kwargs['landcolor'] = 'none'
    kwargs['coast'] = coastline
    kwargs['decimal_degrees'] = True
    fig, ax = cplt.create(extent, **kwargs)

    pargs = dict()
    pargs['norm_clevs'] = norm
    pargs['extend'] = 'neither'
    pargs['cmap'] = cmap
    pargs['clab'] = color_label
    plot_pcolormesh(fig, ax, lon, lat, vsub.values, **pargs)

    if h:
        title = f'{vsub.long_name} {h}m: Wind Farm\nSummer 2022'
        save_file = f'{v}_{h}m_avg-1km_wf2km_nyb-surfacemap.png'
    else:
        title = f'{vsub.long_name}: Wind Farm\nSummer 2022'
        save_file = f'{v}_avg-1km_wf2km_nyb-surfacemap.png'
    #ax.set_title(title, pad=8)

    # shapefiles downloaded from the BOEM website: https://www.boem.gov/renewable-energy/mapping-and-data/renewable-energy-gis-data
    lease = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/BOEM-Renewable-Energy-Shapefiles_11_2_2022/Wind_Lease_Outlines_11_2_2022.shp'
    #lease = glob.glob('/home/coolgroup/bpu/mapdata/shapefiles/BOEM-Renewable-Energy-Shapefiles-current/Wind_Lease_Outlines*.shp')[0]
    oargs = dict()
    oargs['edgecolor'] = 'magenta'
    map_add_boem_outlines(ax, lease, **oargs)

    plt.savefig(os.path.join(save_dir, save_file), dpi=200)

######################################################################################################################
    # plot average: control
    plot_pcolormesh(fig, ax, lon, lat, vctrlsub.values, **pargs)
    map_add_boem_outlines(ax, lease, **oargs)

    if h:
        title = f'{vctrlsub.long_name} {h}m: Control\nSummer 2022'
        save_file = f'{v}_{h}m_avg-1km_ctrl-surfacemap.png'
    else:
        title = f'{vctrlsub.long_name}: Control\nSummer 2022'
        save_file = f'{v}_avg-1km_ctrl-surfacemap.png'

    #ax.set_title(title, pad=8)
    plt.savefig(os.path.join(save_dir, save_file), dpi=200)
    plt.close()

######################################################################################################################
    # plot difference
    # define color map and label
    try:
        cmin = vconfigs['difflims'][0]
        cmax = vconfigs['difflims'][1]
        bins = vconfigs['diffbins']
    except KeyError:
        cmax = np.ceil(np.nanmax(abs(diff)))
        cmin = -cmax
        bins = (cmax - cmin) * 10

    try:
        cmap = vconfigs['diffcmap']
    except KeyError:
        cmap = plt.get_cmap('RdBu_r')

    levels = MaxNLocator(nbins=bins).tick_values(cmin, cmax)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    #norm = TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)

    if diff.units == 'm s-1':
        unitslab = '($\\mathrm{{m}}\\ \\mathrm{{s}}^{{-1}}$)'
    elif diff.units == 'kg kg-1':
        unitslab = '($\\mathrm{{kg}}\\ \\mathrm{{kg}}^{{-1}}$)'
    elif diff.units == 'W m-2':
        unitslab = '($\\mathrm{{W}}\\ \\mathrm{{m}}^{{-2}}$)'
    else:
        unitslab = f'({diff.units})'

    if diff.long_name == 'Average UPWARD HEAT FLUX AT THE SURFACE Difference':
        ln = 'Average Upward Heat Flux Difference'
    else:
        ln = diff.long_name
    color_label = f'{ln} {unitslab}'

    kwargs = dict()
    kwargs['figsize'] = (8, 8)
    kwargs['oceancolor'] = 'none'
    kwargs['landcolor'] = 'white'
    kwargs['coast'] = coastline
    kwargs['decimal_degrees'] = True
    fig, ax = cplt.create(extent, **kwargs)

    pargs = dict()
    pargs['norm_clevs'] = norm
    pargs['extend'] = 'neither'
    pargs['cmap'] = cmap
    pargs['clab'] = color_label

    if v in ['ws', 'ws10']:
        contour_list = [-0.5]
        add_contours(ax, lon, lat, masked_diff, contour_list, label_format='%.1f', color='magenta', linewidth=1, linestyle='--')
    if v in ['UST']:
        pargs['cbar_ticks'] = [-0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04]
        contour_list = [-0.02]
        add_contours(ax, lon, lat, masked_diff, contour_list, label_format='%.2f', color='magenta', linewidth=1, linestyle='--')
        contour_list = [-0.024]
        add_contours(ax, lon, lat, masked_diff, contour_list, label_format='%.2f', color='blue', linewidth=1,
                        linestyle='--')
    # if v in ['TKE_PBL']:
    #     contour_list = [-0.05]
    #     add_contours(ax, lon, lat, masked_diff, contour_list, label_format='%.2f', color='red', linewidth=1, linestyle='--')
    if v in ['TEMP']:
        contour_list = [-.1, .1]
        add_contours(ax, lon, lat, masked_diff, contour_list, label_format='%.1f', color='blue', linewidth=1,
                        linestyle='--')
    if v in ['T2', 'TSK']:
        contour_list = [-.05, .05]
        add_contours(ax, lon, lat, masked_diff, contour_list, label_format='%.2f', color='magenta', linewidth=1,
                        linestyle='--')
    if v in ['HFX']:
        pargs['cbar_ticks'] = [-1.5, -1, -0.5, 0, 0.5, 1, 1.5]
    #     contour_list = [.5]
    #     add_contours(ax, lon, lat, masked_diff, contour_list, label_format='%.2f', color='blue', linewidth=1,
    #                     linestyle='--')
    if v in ['Q2']:
        pargs['cbar_ticks'] = [-0.00015, -0.0001, -0.00005, 0, 0.00005, 0.0001, 0.00015]

    plot_pcolormesh(fig, ax, lon, lat, masked_diff, **pargs)

    oargs = dict()
    oargs['edgecolor'] = 'gray'
    map_add_boem_outlines(ax, lease, **oargs)

    if h:
        title = f'{diff.long_name} {h}m (wind farm - control)\nSummer 2022'
        save_file = f'{v}_{h}m_avg-diff-surfacemap.png'
    else:
        title = f'{diff.long_name} (wind farm - control)\nSummer 2022'
        save_file = f'{v}_avg-diff-surfacemap.png'

    #ax.set_title(title, pad=8)
    if v in ['Q2']:
        plt.subplots_adjust(right=.82, left=0.08)
    plt.savefig(os.path.join(save_dir, save_file), dpi=200)

    plt.close()

######################################################################################################################
    # plot proportion reduction from the control
    # define color map and label
    prop = diff / vctrlsub
    print(f'Proportion reduction from control range: [{str(np.round(np.nanmin(prop), 3))}, {str(np.round(np.nanmax(prop), 3))}]')
    prop.attrs['long_name'] = f'{vsub.attrs["long_name"]} Diff / Ctrl'

    if v == 'TKE_PBL':
        cmax = .5
        cmin = -.5
        bins = 30
    elif v == 'HFX':
        cmax = 1
        cmin = -1
        bins = 30
    else:
        # round up to the nearest decimal
        x = np.nanmax(abs(prop))
        cmax = ceil(x * 100) / 100.0
        cmin = -cmax
        bins = 30

    cmap = plt.get_cmap('RdBu_r')

    levels = MaxNLocator(nbins=bins).tick_values(cmin, cmax)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    color_label = f'{prop.long_name}'

    kwargs = dict()
    kwargs['figsize'] = (8, 8)
    kwargs['oceancolor'] = 'none'
    kwargs['landcolor'] = 'none'
    kwargs['coast'] = coastline
    kwargs['decimal_degrees'] = True
    fig, ax = cplt.create(extent, **kwargs)

    pargs = dict()
    pargs['norm_clevs'] = norm
    pargs['extend'] = 'neither'
    pargs['cmap'] = cmap
    pargs['clab'] = color_label

    if v in ['ws', 'ws10', 'UST']:
        contour_list = [-0.05]
        add_contours(ax, lon, lat, prop, contour_list, label_format='%.2f', color='red', linewidth=1,
                        linestyle='--')
        contour_list = [-0.2]
        add_contours(ax, lon, lat, prop, contour_list, label_format='%.2f', color='cyan', linewidth=1,
                        linestyle='--')
    if v in ['Q2']:
        contour_list = [-0.005]
        add_contours(ax, lon, lat, prop, contour_list, label_format='%.3f', color='red', linewidth=1,
                        linestyle='--')
    # if v in ['TKE_PBL']:
    #     contour_list = [-0.1]
    #     add_contours(ax, lon, lat, prop, contour_list, label_format='%.2f', color='red', linewidth=1,
    #                     linestyle='--')

    plot_pcolormesh(fig, ax, lon, lat, prop, **pargs)

    oargs = dict()
    oargs['edgecolor'] = 'gray'
    map_add_boem_outlines(ax, lease, **oargs)

    if h:
        title = f'{prop.long_name} {h}m\nSummer 2022'
        save_file = f'{v}_{h}m_diff-proportion-surfacemap.png'
    else:
        title = f'{prop.long_name}\nSummer 2022'
        save_file = f'{v}_diff-proportion-surfacemap.png'

    #ax.set_title(title, pad=8)
    plt.savefig(os.path.join(save_dir, save_file), dpi=200)

    plt.close()


if __name__ == '__main__':
    v = 'Q2'  # ws10 ws UST Q2 T2 TEMP TEMP TKE_PBL TKE_PBL HFX TSK
    height = 2  # 10 160 None 2 2 200 120 20 200 None None
    coast = 'full'  # full low
    savedir = '/Users/garzio/Documents/rucool/Miles/RMI/2024_Ocean_Mixing/draft_plots_v3'
    main(v, height, coast, savedir)
