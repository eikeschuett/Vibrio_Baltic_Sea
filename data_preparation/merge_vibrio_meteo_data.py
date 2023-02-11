# -*- coding: utf-8 -*-
"""
Calculates mean and trend with the fexible moving window for the COSMO REA6
meteorological reanalysis data from the DWD.
"""

import os
import xarray as xr
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

def avg_ws_to_nc():
    """ calculate the average wind speed from U and V components and export the
    results in a new netCDF file"""
    nc_dir = "N:/data/meteo_data/COSMO_REA6/netCDF/"
    fpath = os.path.join(nc_dir, "U_10M.2D.DayMean.nc4")
    with xr.open_dataset(fpath) as xds_u:
        fpath = os.path.join(nc_dir, "V_10M.2D.DayMean.nc4")
        with xr.open_dataset(fpath) as xds_v:
            xds_u["ws_mean"] = np.sqrt(xds_u["10u"]**2 + xds_v["10v"]**2)
            xds_u["ws_mean"].attrs = {'standard_name': 'windspeed',
                                      'long_name': '10 metre mean wind speed',
                                      'units': 'm s**-1',
                                      'cell_methods': 'time: mean. '
                                                      'calc: sqrt(u**2+v**2)'}
            export_dir = os.path.join(nc_dir, "WS_10m.2D.DayMean.nc4")
            xds_u["ws_mean"].to_dataset().to_netcdf(export_dir)

def get_var_at_date(xds, time, tolerance = pd.Timedelta(hours = 12)):
    """returns the variable at a date within a specified tolerance"""
    try:
        return xds.sel(time = time, 
                       method = "nearest", 
                       tolerance = tolerance).values.item()
    except KeyError:
        return np.nan
    
def get_data_of_last_month(xds, var, tup):
    """returns the data of the last 30 days and the sampling date. tup is a 
    tuple containing (date, lat, lon)"""
    date = pd.to_datetime(tup[0])
    daterange = pd.date_range(date-pd.Timedelta(days = 30), periods = 30)
    lat = tup[1]
    lon = tup[2]
    
    xds_sub = xds.sel(lat = lat,
                      lon = lon,
                      time = daterange,
                      method = "nearest")[var].drop_duplicates(dim="time")
    
    return [get_var_at_date(xds_sub, time) for time in daterange]

def linreg(X, Y):
    """ Fits a linear model to X and Y data"""
    
    # first, remove potential NaNs from Y
    X = X[~np.isnan(Y)]
    Y = Y[~np.isnan(Y)]
    
    # if enough data remains, fit the linear model and return the slope
    if len(X) > 1:
        reg = LinearRegression()
        reg.fit(X, Y)
        return reg.coef_[0]
    # otherwise return a NaN
    else:
        return np.nan
    
def merge_vibrio_meteo_data(df, var_dict):
    """
    Extracts environmental data for each sampling and each variable and export 
    the data as csv-files into a folder, which is specified within 
    this function.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containg location and datetime of all samples for which the
        environmental data shall be extracted. Expects columns names "date",
        "east" and "north" (with UTM coordinates) or "date", "lat", "lon" (if 
        crs_proj is WGS84).
    var_dict : dict
        Dict containing information on each variable for which data shall be 
        extracted. Needs to have at least the following keys: "filename" 
        (of the netCDF file), "var" (name of the variable in the netCDF file), 
        "offset" (offset which will be applied to the result. May be useful to
        convert K to °C).

    Returns
    -------
    str
        Returns the name of the variable after the computation is completed.
        Main output is stored as a CSV file to "N:/data/merged_vibrio_env/" 
        (see variable "out_fpath" below).

    """
    # specify the period of studied time lags
    dates = np.arange(-29,1)
    
    # open the netCDF file
    nc_dir = "N:/data/meteo_data/COSMO_REA6/netCDF/"
    with xr.open_dataset(os.path.join(nc_dir, var_dict["filename"])) as xds:
        
        # compile the data of the last 30 before each measurement and add it 
        # to df
        a = [get_data_of_last_month(xds, var_dict["var"], tup) for tup in 
             zip(df["date"], df["lat"], df["lon"])]
        
        # create a DataFrame from the result
        df_temp = pd.DataFrame(a, columns = abs(dates))
        # reorder the columns from 0 to (-)30 to be consistent with other codes
        df_temp = df_temp[df_temp.columns[::-1]]
        # apply an offset, if there is one (to convert from Kelvin to °C)
        df_temp = df_temp + var_dict["offset"]
        
    ### calculate statistics (mean and slope of LinReg)
    # prepare empty lists to store the data
    columns_mean = []
    columns_linreg = []
    means = []
    linreg_results = []
    
    for start in np.flip(abs(dates)):
        # if start > 5:
        #     continue
        for end in np.arange(start,30):
            
            columns_mean.append("{}_mean_{}_{}".format(var_dict["var"], 
                                                       start, 
                                                       end))
            
            columns_linreg.append("{}_linreg_{}_{}".format(var_dict["var"], 
                                                           start, 
                                                           end))
            
            # if start and end are not the same dates
            if start != end:
                # get the slice of the original dataframe 
                df_slice = df_temp.loc[:,start:end]
                means.append(df_slice.mean(axis = 1).values)
                
                # linear regression
                linreg_results.append(df_slice.apply(lambda x: linreg(
                                           X = np.arange(start, 
                                                         end+1).reshape(-1, 1),
                                           Y = x.values), 
                                      axis=1).to_numpy())
    
            else:
                # if start and end are the same, simply copy the data as mean
                means.append(df_temp.loc[:,start].values)
                
                # no linreg can be calculated, so insert NaNs
                linreg_results.append(np.repeat(np.nan,len(df_temp)))
    
    out_fpath = "N:/data/merged_vibrio_env/"
    
    ### Create a dataframe from the mean results, concat it with df and save it
    df_mean = pd.DataFrame(means).T
    df_mean.columns = columns_mean
    df_mean = pd.concat([df, df_mean], axis = 1)
    df_mean.to_csv(os.path.join(out_fpath + "vibrio_mean_{}.csv".format(
                                                            var_dict["var"])),
                   sep = "\t", 
                   # index = False,
                   )
    
    # do the same with the result of the linear regression
    df_linreg = pd.DataFrame(linreg_results).T
    df_linreg.columns = columns_linreg
    df_linreg = pd.concat([df, df_linreg], axis = 1)
    df_linreg.to_csv(os.path.join(out_fpath + 
                                  "vibrio_linreg_{}.csv".format(
                                                          var_dict["var"])),
                   sep = "\t", 
                   # index = False,
                   )
    
    return var_dict["var"]
    
#%%
if __name__ == "__main__":
    from tqdm import tqdm
    
    #calculate the average wind speed from U and V components and export the
    # results in a new netCDF file
    avg_ws_to_nc()
    
    # import the CSV containing all in-situ measurements
    df = pd.read_csv("N:/data/vibrio/cleaned/vibrio_vuln_cleaned_sh_mv.csv",
                     sep = "\t")
    
    # create a dict containing the filenames, variables and offsets
    vars_dicts = [{"filename": "ASWDIFD_S.2D.DayMean.nc4",
                   "var": "ASWDIFD_S",
                   "offset": 0},
                  {"filename": "T_2M.2D.DayMean.nc4",
                   "var": "T_2M_CL",
                   "offset": -273.15}, # convert from Kelvin to °C
                  {"filename": "TOT_PRECIP.2D.DaySum.nc4",
                   "var": "tp",
                   "offset": 0},
                  {"filename": "WS_10m.2D.DayMean.nc4",
                   "var": "ws_mean",
                   "offset": 0},
                   ]
    
    # run the function to extract environmental data for each sampling and each
    # variable and export the data as csv-files into a folder, which is 
    # specified within merge_vibrio_meteo_data()
    [merge_vibrio_meteo_data(df, var_dict) for var_dict in tqdm(vars_dicts)]


