# -*- coding: utf-8 -*-
"""
Class to create match ups between models and in-situ point data. It calculates 
error statistics and provides a plot-method to visualize results in a 
scatterplot. It is applied in `validate_GETM_SST_SSS.py`.
"""

import pandas as pd
import numpy as np
# import rioxarray
# import xarray as xr
import matplotlib.pyplot as plt
# import datetime
from datetime import timedelta
from shapely.geometry import Point
import pyproj# import Proj, transform 
import warnings

# import cartopy.crs as ccrs
# import pyogrio as pg
from math import ceil

def calc_distance_wgs84(p1, p2, source_crs, target_epsg = 32632):
    proj = pyproj.Transformer.from_crs(source_crs, 
                                       target_epsg, 
                                       always_xy=True)
    
    p1 = Point(proj.transform(p1.x, p1.y))
    p2 = Point(proj.transform(p2.x, p2.y))
    
    return p1.distance(p2)

class match_up(object):
    def __init__(self, gdf, xds, time_col, #data_col, 
                 #data_std = None,
                 #xds_var = "data", 
                 var_dict,
                 station = None,
                 alpha = None,
                 markersize = None,
                 xds_offset = 0.0,
                 max_td = timedelta(hours = 24),
                 max_dist = 5000,
                 verbose = 0):
        """
        Initializes a new insatnce of the match_up-class.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            GeoDataFrame containing the in-situ data.
        xds : xarray.Dataset
            Xarray dataset which shall be validated. Needs to have a CRS 
            assigned using rioxarray (e.g. xds.rio.write_crs(4326, 
                                                             inplace=True))
            The xds is expected to have variables called "lat" and "lon".
        time_col : string
            Name of the column in gdf containnig pd.Datatime objects (or 
            similar) indicating the in-situ measurement time.
        var_dict : dict
            Dictionary of the columns in gdf and the variables in xds to be 
            compared. Expects columns names in gdf as keys and xds variable 
            names as values. E.g. {"ctdtmp": "temp", "sal": "salt"}. 
            Standard deviations of the data in gdf can be passed if the column
            name is "variable_std".
        data_std : string, optional
            Name of the column in gdf containing the standard deviation of the
            of the measured parameter. The default is None.
        xds_var : string, optional
            Name of the variable in xds which shall be validated. The default
            is "data".
        xds_offset : float, optional
            Offset to be applied on the xds-data (e.g. to convert Kelvin to 
            °C.) The default is 0.0.
        station: list or column, optional
            List that contains the stations of the in-situ data. Should be the 
            same length as gdf. 
        alpha: list or column, optional
            List that contains the preferred alpha for the scatterplot of the 
            match-up analysis. Should be the same length as gdf.
        markersize: list or column, optional 
            List for markersize of each point in scatterplot. Should be the 
            same length as gdf.
        max_td : datetime.timedelta, optional
            Maximum temporal offset between in-situ measurement and modelled 
            time in xds. The default is timedelta(hours = 24).
        max_dist : int/float, optional
            Maximum tolerated distance between in-situ measurement position and
            pixel coordinate in xds [in meters!]. Should be larger than the 
            spatial resolution of the grid. The default is 5000.
        verbose : float, optional
            Set to 1 to raise a warning if the timedelta or distance thresholds
            are exceeded. The default is 0.
        """
        
        # if necessary, reproject the gdf to the crs of the xds
        if gdf.crs != xds.rio.crs:
            gdf = gdf.to_crs(xds.rio.crs)

        # prepare empty dicts to store results in
        self.data = dict()
        self.rmse = dict()
        self.mae = dict()
        self.bias = dict()
        self.n = dict()
        self.corrcoef = dict()
        
        # iterate over each varibale which shall be compared
        for gdf_col, xds_var in var_dict.items():
            # create the name of a potential std-column
            gdf_col_std = gdf_col + "_std"
            
            # empty list to store all match ups for this variable
            match_ups = []
            
            # iterate over each in-situ observation
            for i, row in gdf.iterrows():
                # find the nearest modelled value to the observation
                xds_point = xds.sel(lat = row.geometry.y, 
                                    lon = row.geometry.x,
                                    time = row[time_col],
                                    method = "nearest")[xds_var]
                
                # check if the time difference between observation and model is
                # within the specified limits
                td = abs(row[time_col] - xds_point.time.values)
                if td <= max_td:
                    
                    # calculate distance between observation and grid cell 
                    xds_pt_geom = Point([xds_point["lon"], 
                                         xds_point["lat"]])
                    
                    # get the distance in meters
                    distance = calc_distance_wgs84(p1 = xds_pt_geom, 
                                                   p2 = row["geometry"],
                                                   source_crs = xds.rio.crs,
                                                   target_epsg = 32632)
                    if (markersize is None):
                        row[markersize] = 10
                    
                    if (alpha is None):
                        row[alpha] = 1
                        
                    if (station is None):
                        row[station] = "black"
                    
                    # if the distance is within the limits, add the match up
                    if distance <= max_dist:
                        # preserve information on std of observed value if it 
                        # exists
                        if gdf_col_std in gdf.columns:
                            match_ups.append({
                                "time_obs": row[time_col],
                                "latitude": row.geometry.y,
                                "longitude": row.geometry.x,
                                "time_model": xds_point.time.values,
                                "distance": distance,
                                "modelled": float(xds_point.values) + 
                                            xds_offset,
                                "observed": row[gdf_col],
                                "residual": float(xds_point.values) + 
                                            xds_offset - 
                                            row[gdf_col],
                                "observed_std": row[gdf_col_std]
                                            })
                        else:
                            match_ups.append({
                                "time_obs": row[time_col],
                                "latitude": row.geometry.y,
                                "longitude": row.geometry.x,
                                "time_model": xds_point.time.values,
                                "distance": distance,
                                "modelled": float(xds_point.values) + 
                                            xds_offset,
                                "observed": row[gdf_col],
                                "residual": float(xds_point.values) + 
                                            xds_offset - 
                                            row[gdf_col],
                                "station":  row[station],
                                "alpha":  row[alpha],
                                "markersize":  row[markersize]
                                                              })
                       
                                
                    else:
                        if verbose == 1:
                            warnings.warn("Closest match-up with row {} "
                                          "exceeds the max_dist limit. "
                                          "Distance between Point and pixel is"
                                          " {}.".format(i, distance))
                        
                else:
                    if verbose == 1:
                        
                        warnings.warn("Closest match-up with row {} exceeds "
                                      "the maximum time difference limit. "
                                      "Temporal offset between both "
                                      "measurements is {}.".format(i, td))
            
            # convert match ups to a dataframe
            df_tmp = pd.DataFrame(match_ups)
            
            self.data[xds_var] = df_tmp
            
            # calculate statistics
            self.rmse[xds_var] = ((df_tmp["modelled"] - 
                                   df_tmp["observed"]) ** 2).mean() ** 0.5
            self.mae[xds_var] = (abs(df_tmp["modelled"] - 
                                     df_tmp["observed"])).mean()
            self.bias[xds_var] = (df_tmp["modelled"] - 
                                  df_tmp["observed"]).mean()
            self.n[xds_var] = df_tmp["residual"].count()
            self.corrcoef[xds_var] = df_tmp["modelled"].corr(
                                            df_tmp["observed"])
    
    #Plot function 
    #! Important: markersize, alpha and colour needs to be set in the match-up 
    # function! (Variables:"station","alpha" and "markersize")
    def plot(self, xlabel = "observed", ylabel = "modelled",
             nrows = "auto", ncols = 3, figsize = (5,5),
             ax = None, key = None,
             title = None, 
             plot_stats = False,
             plot_type = "scatterplot",
             
             save_as = None):
        
        
        
        if nrows == "auto":
            if ncols > len(self.data):
                ncols = len(self.data)
                nrows = 1
            else:
                nrows = ceil(len(self.data)/ncols)
        
        if not ax and not key:
            
            fig = plt.figure(figsize = figsize)
            
            if title:
                fig.suptitle(title)
            
        user_ax = ax
        
        if all(v.empty for (k,v) in self.data.items()):
            raise IndexError("No valid match-ups that can be plotted.")
            
        else:
            
            for i, (k, v) in enumerate(self.data.items()):
            
                if not user_ax:
                    ax = fig.add_subplot(nrows, ncols, i+1)    
                    if len(self.data.items()) > 1:
                        ax.set_title(k)

                else:
                    if not key == k:
                        continue
            
                ax.axline((0, 0), 
                                slope=1, 
                                color="black", 
                                linestyle="--")
                
                if not v.empty:
                    
                    if plot_type in ["scatter", "scatterplot"]:
                        ax.scatter(v["observed"], v["modelled"], 
                                       alpha = v["alpha"], c = v["station"], 
                                       s = v["markersize"])
                    elif plot_type in ["hexbin", "hexbin_plot"]:
                        hb = ax.hexbin(v["observed"], v["modelled"], 
                                       gridsize=50, cmap='inferno')
                        fig.colorbar(hb, ax = ax, label='N')
                    else:
                        raise ValueError("plot_type must be one of 'scatter' "
                                         "or 'hexbin'.")
                    if "observed_std" in v.columns:
                        ax.errorbar(v["observed"], v["modelled"], 
                                     yerr = v["observed_std"], fmt = "none")
        
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.set_aspect("equal")
                max_lim = np.max([ax.get_ylim()[1], ax.get_xlim()[1]])
                ax.set_ylim([0, max_lim])
                ax.set_xlim([0, max_lim])
                
                if plot_stats:
                    ax.text(0.05, 0.82, '$r$: {:.2f}\n'
                                        #'RMSE: {:.2f}\n'
                                        'MAE: {:.2f}\n'
                                        'bias: {:.2f}\n'
                                        'n: {}'.format(self.corrcoef[k],
                                                       #self.rmse[k],
                                                       self.mae[k],
                                                       self.bias[k],
                                                       self.n[k]).replace("-", 
                                                                  u"\u2212"), 
                            horizontalalignment='left',
                            verticalalignment='center', 
                            transform=ax.transAxes)
            
        if save_as:
            plt.savefig(save_as, dpi = 150, bbox_inches = "tight")
        
        if not user_ax:
            plt.show()
            
            


