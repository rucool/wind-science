#! /usr/bin/env python

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from mpl_toolkits.axes_grid1 import make_axes_locatable


def map_add_boem_outlines(ax, shpfile, edgecolor=None, projection=None, zorder=None, alpha=None):
    edgecolor = edgecolor or 'black'
    projection = projection or ccrs.PlateCarree()
    zorder = zorder or 1
    alpha = alpha or 1

    shape_feature = cfeature.ShapelyFeature(
        Reader(shpfile).geometries(),
        projection,
        edgecolor=edgecolor,
        facecolor='none'
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
