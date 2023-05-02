#! /usr/bin/env python

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from mpl_toolkits.axes_grid1 import make_axes_locatable


def add_contours(ax, londata, latdata, vardata, clist, label_format=None):
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

    CS = ax.contour(londata, latdata, vardata, clist, colors='black', linewidths=.5, transform=ccrs.PlateCarree())
    ax.clabel(CS, inline=True, fontsize=10.5, fmt=label_format)


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


def plot_contourf(fig, ax, x, y, c, cmap, levels=None, ttl=None, clab=None, cbar_ticks=None, extend='both',
                  shift_subplot_right=0.88, xlab=None, ylab=None, yticks=None):
    """
    Create a filled contour plot with user-defined levels and colors
    :param fig: figure object
    :param ax: plotting axis object
    :param x: x-axis data
    :param y: y-axis data
    :param c: color data
    :param cmap: colormap
    :param levels: optional list of data levels
    :param ttl: optional plot title
    :param clab: optional colorbar label
    :param cbar_ticks: optional, specify colorbar ticks
    :param extend: optional, different colorbar extensions, default is 'both'
    :param shift_subplot_right: optional, specify shifting the subplot, default is 0.88
    :param xlab: optional x-label
    :param ylab: optional y-label
    :param yticks: specify optional yticks
    :returns fig, ax objects
    """

    plt.subplots_adjust(right=shift_subplot_right)
    if ttl:
        plt.title(ttl, fontsize=17)
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    fig.add_axes(cax)

    if levels:
        try:
            cs = ax.contourf(x, y, c, levels=levels, cmap=cmap, extend=extend, transform=ccrs.PlateCarree())
        except ValueError:
            cs = ax.contourf(x, y, c, levels=levels, cmap=cmap, extend=extend)
    else:
        try:
            cs = ax.contourf(x, y, c, cmap=cmap, extend=extend, transform=ccrs.PlateCarree())
        except ValueError:
            cs = ax.contourf(x, y, c, cmap=cmap, extend=extend)

    if cbar_ticks:
        cb = plt.colorbar(cs, cax=cax, ticks=cbar_ticks)
    else:
        cb = plt.colorbar(cs, cax=cax)

    if clab:
        cb.set_label(label=clab, fontsize=14)

    if xlab:
        ax.set_xlabel(xlab)
    if ylab:
        ax.set_ylabel(ylab)
    if yticks:
        ax.set_yticks(yticks)

    return fig, ax


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
    if cbar_ticks != 'none':
        cb.ticks = cbar_ticks


def plot_regions():
    regions = dict(
        full_grid=dict(quiver_subset=dict(_10m=11, _160m=13, _200m=13, _250m=13, _500m=13),
                       quiver_scale=45,
                       extent=[-79.79, -69.2, 34.5, 43],
                       xticks=[-78, -76, -74, -72, -70],
                       yticks=[36, 38, 40, 42],
                       subset=False,
                       lease_area=False),
        mab=dict(quiver_subset=dict(_10m=7, _160m=8, _200m=8, _250m=8, _500m=8),
                 quiver_scale=40,
                 extent=[-77.2, -69.6, 36, 41.8],
                 xticks=[-75, -73, -71],
                 yticks=[37, 39, 41],
                 subset=True,
                 lease_area=True),
        nj=dict(quiver_subset=dict(_10m=4, _160m=5, _200m=5, _250m=5, _500m=5),
                quiver_scale=40,
                extent=[-75.7, -71.5, 38.1, 41.2],
                xticks=[-75, -74, -73, -72],
                yticks=[39, 40, 41],
                subset=True,
                lease_area=True),
        southern_nj=dict(quiver_subset=dict(_10m=3, _160m=3, _200m=3, _250m=3, _500m=3),
                         quiver_scale=40,
                         extent=[-75.6, -73, 38.6, 40.5],
                         xticks=[-75, -74.5, -74, -73.5],
                         yticks=[39, 39.5, 40],
                         subset=True,
                         lease_area=True),
        windturb=dict(quiver_subset=dict(_10m=3, _160m=3, _200m=3, _250m=3, _500m=3),
                      quiver_scale=40,
                      extent=[-74.8, -73.7, 38.8, 39.8],
                      xticks=[-74.7, -74.4, -74.1, -73.8],
                      yticks=[39, 39.3, 39.6],
                      subset=True,
                      lease_area=True)
    )

    return regions
