# -*- coding: utf-8 -*-
"""
Here, a function for plotting maps of the study area is defined as well as a
function to remove the Schlei and Bodden areas, where the GETM SSS product
showed poor resulst. Moreover, the map showing the study area with Vibrio 
sampling locations (Figure 1) and the climatologies of the GETM output 
(Figure 4) are created.
"""

import geopandas as gpd    
import xarray as xr
import rioxarray
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import AnchoredText

# 
def remove_Schlei(xds, xdim = "lon", ydim= "lat", 
                  fpath = "N:/data/bbox/Exclude_zone.geojson"):
    """
    helper function to remove Schlei and Bodden from dataset. Clips an 
    xarray.dataset to the extent of a geojson stored in fpath.
    """

    gdf = gpd.read_file(fpath)

    xds = xds.rio.set_spatial_dims(xdim, ydim, inplace=True)
    clipped = xds.rio.clip(gdf.geometry.values, ccrs.UTM(zone = 32),
                           invert=True)
    return clipped
#

def create_map(fig = None, sp_x = 1, sp_y = 1, sp_n = 1,
               gridticksize = False,
               left_labels = True,
               bottom_labels = True,
               extent = [525000, 850000, 5940000, 6086000],
               figsize = (12,10)):
    """
    Creates a geoaxis (or a Figure, if none is provided) of the study region,
    which can be used to easily plot additional data to it.
    """

    # use UTM projection 32 for the plot
    utm_crs = ccrs.UTM(zone = 32)
    
    # If the figure already exists, add a suplot
    if fig:
        ax = fig.add_subplot(sp_x,sp_y,sp_n, projection=utm_crs)
    else:
        # otherwise create the figure first
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(sp_x,sp_y,sp_n, projection=utm_crs)

    # set the extent to the study region
    ax.set_extent(extent,
                  crs=utm_crs)
    
    # add land as background in gray
    gdf = gpd.read_file("N:/data/bbox/EU_Boarders_NUTS_0.geojson")
    gdf = gdf[gdf["CNTR_CODE"].isin(["DE", "DK", "PL"])]
    
    for geom in gdf["geometry"]:
        shape_feature = cfeature.ShapelyFeature([geom], 
                                                ccrs.PlateCarree(), 
                                                facecolor="gray", 
                                                alpha = 0.5,
                                                edgecolor='black', 
                                                lw=1,
                                                zorder = -50)
        ax.add_feature(shape_feature)
    
    # add grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.ylocator = mticker.FixedLocator([54, 54.5])
    gl.xlocator = mticker.FixedLocator([10, 12, 14])
    gl.left_labels = left_labels
    gl.bottom_labels = bottom_labels
    if gridticksize:
        gl.xlabel_style = {'size': gridticksize}
        gl.ylabel_style = {'size': gridticksize}
    
    return fig, ax

def add_subplot_char(ax, label, loc = "upper right"):
    """
    Helper function to draw subplot labels into the corners of the plot.
    """
    text  = AnchoredText("({})".format(label), 
                  loc = loc, 
                  frameon=True)
    
    ax.add_artist(text)

def add_aligned_colorbar(fig, ax, pc, label = None, 
                         ticks = None, ticklabels = None,
                         size = "3%", pad = 0.1):
    """
    Helper function to create a colorbar which is aligned to the axis of a 
    cartopy map.
    """
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="3%", pad=0.1, axes_class=plt.Axes)
    
    fig.add_axes(ax_cb)
    cbar = plt.colorbar(pc, cax=ax_cb, label = label, ticks = ticks)
    
    if type(ticklabels) == list:
        cbar.ax.set_yticklabels(ticklabels)
    
    return cbar, ax_cb


  
  #%% prepare the data
if __name__ == "__main__":
    
    import os
    os.chdir("N:/data/")
    import pandas as pd
    import numpy as np
    import cartopy
    from cmocean import cm
    import cartopy.mpl.geoaxes
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from shapely import geometry
    import matplotlib.patches as mpatches
    
    #%% plot of the study region with the Vibrio samples (Figure 1)
    df = pd.read_csv("merged_vibrio_env/vibrio_data_merged_sh_mv.csv",
                     sep = "\t")
    
    df["month"] = [x.split("-")[1] for x in df.date]
    df["quantity"][df["month"] != "08"] = np.nan
    
    df = df[["station_id", "station_name", "quantity", "lat", "lon"]]
    
    
    df = df.groupby('station_id') \
                       .agg({'station_id':'size', 'quantity':'mean', 
                             'lat': 'mean', 'lon': 'mean'}) \
                       .rename(columns={'station_id':'counts',
                                        'quantity':'quantity_mean'}) \
                       .reset_index()
    
    #% plot the data
    fig, ax = create_map(sp_x = 1, 
                         sp_y = 1, 
                         sp_n = 1,
                         extent = [525000, 850000, 5940000, 6099000],
                         figsize = (12,10))
    
    pc = plt.scatter(x=df.lon, y=df.lat,
                     c = df.counts,
                     cmap = cm.dense,
                     vmin = 1,
                     vmax = 100,
                     alpha = 1,
                     edgecolor = "k",
                     label='_nolegend_',
                     transform = ccrs.PlateCarree())
       
    be = plt.scatter(y=54.5295, x=10.0393,
                     c = "tab:blue",
                     cmap = cm.dense,
                     marker = "P",
                     s = 100,
                     alpha = 1,
                     label = "Boknis Eck Station",
                     edgecolor = "k",
                     transform = ccrs.PlateCarree())
    
    be = plt.scatter(y=54.156053, x=11.763174,
                     c = "tab:green",
                     cmap = cm.dense,
                     marker = "P",
                     s = 100,
                     alpha = 1,
                     label = "Kühlunsgborn Station",
                     edgecolor = "k",
                     transform = ccrs.PlateCarree())
        
    ab = plt.scatter(y=54.8833, x=13.8667,
                     c = "tab:orange",
                     cmap = cm.dense,
                     marker = "P",
                     s = 100,
                     alpha = 1,
                     label = "Arkona Basin Station",
                     edgecolor = "k",
                     transform = ccrs.PlateCarree())    
    
    plt.legend(loc="lower right")
    
    cbar, ax_cb = add_aligned_colorbar(fig = fig, ax = ax, pc = pc, 
                                       label = "Vibrio measurements [counts]", 
                                       # ticks = ticks_needed,
                                       # ticklabels = ticks_wanted,
                                       size = "3%", 
                                       pad = 0.1)
   
    # add inset with map of Germany
    axins = inset_axes(ax, width="55%", height="52%", loc="lower left", 
                   bbox_to_anchor=(-0.18,-0.0,1,1), 
                   bbox_transform=ax.transAxes,
                   axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                   axes_kwargs=dict(map_projection=ccrs.UTM(zone = 32)))
    axins.add_feature(cartopy.feature.COASTLINE)
    axins.add_feature(cfeature.BORDERS)
    
    # add land as background in gray
    gdf = gpd.read_file("N:/data/bbox/EU_Boarders_NUTS_0.geojson")
    gdf = gdf[gdf["CNTR_CODE"].isin(["DE", "DK", "PL", "SE", "NL", "BE", "FR", 
                                     "LU", "CH", "CZ", "AT"])]
    
    for geom in gdf["geometry"]:
        shape_feature = cfeature.ShapelyFeature([geom], 
                                                ccrs.PlateCarree(), 
                                                facecolor="gray", 
                                                alpha = 0.5,
                                                edgecolor='black', 
                                                lw=1,
                                                zorder = -50)
    
        axins.add_feature(shape_feature)
    
    axins.set_extent([5.5,15.5,47.2,56],
                  ccrs.PlateCarree()) #For GER centered Minimap
    
    geom = geometry.box(9.38027949524113,
                        53.760717846561995, 
                        14.446971466429765,  
                        54.920370263065934)
    
    axins.add_geometries([geom], ccrs.PlateCarree(), facecolor='none',
                              edgecolor='tab:red', linewidth=2)
    
    plt.savefig("N:/plots/map_vibrio_samples.jpeg",
                dpi = 300,
                bbox_inches = "tight")
    
    plt.show()
    
    #%% plot climatologies of the output of GETM
    
    xds_fpath = "N:/data/GETM/"   
    
    sst = xr.open_dataset(xds_fpath+"mean_sst_jja.nc").rio.write_crs(
                        ccrs.UTM(zone = 32), inplace=True)["temp"]
    sal = xr.open_dataset(xds_fpath+"mean_salinity_2008-2021.nc"
                          
                          ).rio.write_crs(ccrs.UTM(zone = 32), 
                                          inplace=True)["salt"]

    extent = [525000, 850000, 5966530, 6086000]
    figsize = (8,6)
    fig, ax = create_map(sp_x = 2, 
                         sp_y = 1, 
                         sp_n = 1, 
                         figsize = figsize,
                         extent = extent,
                         bottom_labels = False,
                         )
    
    mc = sst.plot(ax=ax,
                        transform=ccrs.PlateCarree(),
                        vmin = sst.min(),
                        vmax = sst.max(),
                        cmap = cm.thermal,
                        add_colorbar = False,
                        label='_nolegend_',
                        zorder = -25)
    
    plt.title(" ")
    
    cbar1, ax_cb1 = add_aligned_colorbar(fig = fig, ax = ax, pc = mc, 
                                        label = "$\overline{SST}_{JJA}$ [°C]" , 
                                        size = "3%", 
                                        pad = 0.1)
    
    excl = gpd.read_file("N:/data/bbox/Exclude_zone.geojson").geometry
    ax.add_geometries(excl, 
                      ccrs.PlateCarree(), 
                      facecolor='none',
                      edgecolor='k',
                      hatch='\\\\\\', 
                      linewidth=1, 
                      alpha = 0.75)
    
    proxy_artist = mpatches.Rectangle((0, 0), 1, 0.1, 
                                      linewidth=1,
                                      hatch='\\\\\\', 
                                      edgecolor='k',
                                      alpha = 0.75,
                                      facecolor='none')

    # manually add the labels here
    ax.legend([proxy_artist], 
              ['Excluded zones'], 
              loc='lower left', 
              fancybox=True)
    
    add_subplot_char(ax, "a")
    
    _, ax = create_map(fig = fig,
                       sp_x = 2, 
                       sp_y = 1, 
                       sp_n = 2, 
                       figsize = figsize,
                       extent = extent,
                       )
    
    mc = sal.plot(ax=ax,
                        transform=ccrs.PlateCarree(),
                        vmin = sal.min(),
                        vmax = sal.max(),
                        cmap = cm.haline,
                        add_colorbar = False,
                        label='_nolegend_',
                        zorder = -25)
    
    plt.title(" ")
    
    cbar1, ax_cb1 = add_aligned_colorbar(fig = fig, ax = ax, pc = mc, 
                                        label = "$\overline{SSS}$ [PSU]" , 
                                        size = "3%", 
                                        pad = 0.1)
    
    geom = geometry.box(9.38027949524113,
                        53.760717846561995, 
                        14.446971466429765,  
                        54.920370263065934)
    axins.add_geometries([geom], 
                         ccrs.PlateCarree(), 
                         facecolor='none',
                         edgecolor='tab:red', 
                         linewidth=2)
    
    excl = gpd.read_file("N:/data/bbox/Exclude_zone.geojson").geometry
    ax.add_geometries(excl, 
                      ccrs.PlateCarree(), 
                      facecolor='none',
                      edgecolor='k',
                      hatch='\\\\\\', 
                      linewidth=1, 
                      alpha = 0.75)
    
    
    
    proxy_artist = mpatches.Rectangle((0, 0), 1, 0.1, 
                                      linewidth=1,
                                      hatch='\\\\\\', 
                                      edgecolor='k',
                                      alpha = 0.75, 
                                      facecolor='none')

    # And manually add the labels here
    ax.legend([proxy_artist], 
              ['Excluded zones'], 
              loc='lower left', 
              fancybox=True)
    
    add_subplot_char(ax, "b")
    
    
    plt.subplots_adjust(hspace = 0.01)
    
    
    plt.savefig("N:/plots/GETM_sal_sst_mean.jpeg",
                dpi = 300,
                bbox_inches = "tight")
    
    plt.show()
