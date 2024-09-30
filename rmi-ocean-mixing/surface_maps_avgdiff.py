#!/usr/bin/env python

"""
Author: Lori Garzio on 9/26/2024
Last modified: 9/26/2024

"""

import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import cmocean as cmo
import cool_maps.plot as cplt
import functions.plotting as pf
plt.rcParams.update({'font.size': 12})  # all font sizes are 12 unless otherwise specified


def subset_grid(ext, data, lon_name, lat_name):
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


def main(variable, h, coastline, plot_region, save_dir):
    plt_region = pf.plot_regions()
    extent = plt_region[plot_region]['extent']

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
        ws=dict(cmap=plt.get_cmap('BuPu'), diffmask=[-0.1, 0.1]),
        ws10=dict(cmap=plt.get_cmap('BuPu'), diffmask=[-0.1, 0.1]),
        UST=dict(cmap=plt.get_cmap('YlGnBu'), diffmask=[-0.01, 0.01]),
        Q2=dict(cmap=cmo.cm.rain, diffmask=[-0.001, 0.001]),
        T2=dict(cmap=cmo.cm.thermal, diffmask=[-0.01, 0.01]),
        TEMP=dict(cmap=cmo.cm.thermal, diffmask=[-0.01, 0.01]),
        TKE_PBL=dict(cmap=cmo.cm.tempo, diffmask=[-0.01, 0.01]),
        HFX=dict(cmap=cmo.cm.thermal, diffmask=[-0.1, 0.1])
    )

    vconfigs = configs[v]
    vsub, lon, lat = subset_grid(extent, ds[varnameavg], 'XLONG', 'XLAT')
    vctrlsub, lon, lat = subset_grid(extent, dsc[varnameavg], 'XLONG', 'XLAT')

    if v in ['TEMP', 'T2']:
        ln = vsub.attrs['long_name']
        vsub = vsub - 273.15  # convert from K to degrees C
        vsub.attrs['units'] = 'deg_C'
        vsub.attrs['long_name'] = ln
        vctrlsub = vctrlsub - 273.15  # convert from K to degrees C
        vctrlsub.attrs['units'] = 'deg_C'
        vctrlsub.attrs['long_name'] = ln

    diff = vsub - vctrlsub
    diff.attrs['units'] = vsub.attrs['units']
    diff.attrs['long_name'] = f'{vsub.attrs["long_name"]} Difference'
    masked_diff = np.ma.masked_inside(diff, vconfigs['diffmask'][0], vconfigs['diffmask'][1])

    # define color map and label
    if v in ['UST']:
        cmin = 0.1
        cmax = 0.3
        bins = 20
    elif v in ['Q2']:
        cmin = 0
        cmax = 0.02
        bins = 20
    elif v in ['TKE_PBL']:
        cmin = 0
        cmax = 2
        bins = 20
    elif v in ['HFX']:
        cmin = -10
        cmax = 10
        bins = (cmax - cmin) * 5
    else:
        cmin = np.floor(np.nanmin([vsub, vctrlsub]))
        cmax = np.ceil(np.nanmax([vsub, vctrlsub]))
        bins = (cmax - cmin) * 5
    cmap = vconfigs['cmap']
    levels = MaxNLocator(nbins=bins).tick_values(cmin, cmax)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    color_label = f'{vsub.long_name} ({vsub.units})'

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
    pf.plot_pcolormesh(fig, ax, lon, lat, vsub.values, **pargs)

    if h:
        title = f'{vsub.long_name} {h}m: Wind Farm\nSummer 2022'
        save_file = f'{v}_{h}m_avg-1km_wf2km_nyb-surfacemap.png'
    else:
        title = f'{vsub.long_name}: Wind Farm\nSummer 2022'
        save_file = f'{v}_avg-1km_wf2km_nyb-surfacemap.png'
    ax.set_title(title, pad=8)

    lease = '/Users/garzio/Documents/rucool/bpu/wrf/lease_areas/BOEM-Renewable-Energy-Shapefiles_11_2_2022/Wind_Lease_Outlines_11_2_2022.shp'
    #lease = glob.glob('/home/coolgroup/bpu/mapdata/shapefiles/BOEM-Renewable-Energy-Shapefiles-current/Wind_Lease_Outlines*.shp')[0]
    oargs = dict()
    oargs['edgecolor'] = 'magenta'
    pf.map_add_boem_outlines(ax, lease, **oargs)

    plt.savefig(os.path.join(save_dir, save_file), dpi=200)

    # plot average: control
    pf.plot_pcolormesh(fig, ax, lon, lat, vctrlsub.values, **pargs)
    pf.map_add_boem_outlines(ax, lease, **oargs)

    if h:
        title = f'{vctrlsub.long_name} {h}m: Control\nSummer 2022'
        save_file = f'{v}_{h}m_avg-1km_ctrl-surfacemap.png'
    else:
        title = f'{vctrlsub.long_name}: Control\nSummer 2022'
        save_file = f'{v}_avg-1km_ctrl-surfacemap.png'

    ax.set_title(title, pad=8)
    plt.savefig(os.path.join(save_dir, save_file), dpi=200)
    plt.close()

    # plot difference
    # define color map and label
    if v in ['UST']:
        cmax = .1
        cmin = -cmax
        bins = 20
    elif v in ['Q2']:
        cmax = .05
        cmin = -cmax
        bins = 20
    elif v in ['T2', 'TEMP']:
        cmax = .2
        cmin = -cmax
        bins = 20
    else:
        cmax = np.ceil(np.nanmax(abs(diff)))
        cmin = -cmax
        bins = (cmax - cmin) * 10

    cmap = plt.get_cmap('RdBu_r')
    levels = MaxNLocator(nbins=bins).tick_values(cmin, cmax)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    color_label = f'{diff.long_name} ({diff.units})'

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

    if v in ['ws', 'ws10']:
        contour_list = [-0.5]
        pf.add_contours(ax, lon, lat, masked_diff, contour_list, label_format='%.1f', color='red', linewidth=1, linestyle='--')

    pf.plot_pcolormesh(fig, ax, lon, lat, masked_diff, **pargs)

    oargs = dict()
    oargs['edgecolor'] = 'gray'
    pf.map_add_boem_outlines(ax, lease, **oargs)

    if h:
        title = f'{diff.long_name} {h}m (wind farm - control)\nSummer 2022'
        save_file = f'{v}_{h}m_avg-diff-surfacemap.png'
    else:
        title = f'{diff.long_name} (wind farm - control)\nSummer 2022'
        save_file = f'{v}_avg-diff-surfacemap.png'

    ax.set_title(title, pad=8)
    plt.savefig(os.path.join(save_dir, save_file), dpi=200)

    plt.close()


if __name__ == '__main__':
    v = 'HFX'  # ws10 ws UST Q2 T2 TEMP TEMP TKE_PBL TKE_PBL HFX
    height = 20  # 10 160 None 2 2 200 120 20 200 20
    coast = 'full'  # full low
    region = 'windturb_nyb'
    savedir = '/Users/garzio/Documents/rucool/Miles/RMI/2024_Ocean_Mixing/draft_plots_v1'
    main(v, height, coast, region, savedir)
