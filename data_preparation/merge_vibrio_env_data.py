# -*- coding: utf-8 -*-
"""
Merges Vibrio data with all biogeochemical and hydrodynamic data containing all
time windows (i.e. up to 30 days before Vibrio measurement).
Note: This code is rather inefficient. It works perfectly, but it may be worth 
to revise it if you need to run it again...
"""

import pandas as pd
import xarray as xr
import numpy as np
from pyproj import CRS, Transformer
from tqdm import tqdm
import warnings
from sklearn.linear_model import LinearRegression

def calc_env_data_offset(fpath, fnames, variable, df, column, days_diff, 
                         function, crs_4326, crs_proj, offset = 0, 
                         init_dist = 0, max_dist = 0, increment = 0):
    """
    Appends environmental data for each moving window of flexible size to the
    dataframe 'df'.

    Parameters
    ----------
    fpath : str
        Path to the location of the netCDF-file.
    fnames : str/list
        String or list containing the filename(s) of the netCDFs from which the
        data shall be extracted.
    variable : str
        Name of the variable in the netCDF file which will be extracted.
    df : pd.DataFrame
        DataFrame containg location and datetime of all samples for which the
        environmental data shall be extracted. Expects columns names "date",
        "east" and "north" (with UTM coordinates) or "date", "lat", "lon" (if 
        crs_proj is WGS84).
    column : str
        Name of the column which will be added to df and populated with the 
        data.
    days_diff : List
        List of date offsets for which the data will be calculated. E.g. if 
        np.arange(0,10,2) is provided, the environmental data of the sample 
        date (0), days 0-1, 0-2, ..., 0-10, 1-1, 1-2, ..., 1-10, 2-2, 2-3, ...,
        2-10, ..., 10-10 will be calculated.
    function : str
        Function by which the data will be aggregated. Must be one of ["mean",
        "min", "max", "median", "sum", "linreg"]. If "linreg" is chosen, the
        slope of a linear regression performed with sklearns LinearRegression
        is returend.
    crs_4326 : CRS
        CRS of the netCDF file.
    crs_proj : CRS
        CRS of the data in df.
    offset : float, optional
        Offset which will be added to the results (useful to recalculate Kelvin
        to °C). The default is 0.
    init_dist : int, optional
        Initial search distance for no-nan values in the raster. The 
        unit is the same as the coordinates in df.east and df.north (i.e. m if 
        coordinates are in a metric system. The default is 0.
    max_dist : int, optional
        Maximum search distance for no-NaN values in the raster. Same unit as 
        init_dist. The default is 0.
    increment : float, optional
        Increment by which the search window is widened after each iterarion. 
        Same unit as init_dist. The default is 0.

    Returns
    -------
    df : pd.DataFrame
        DataFrame containing the initial data of df. Environmental data is 
        appended to new columns.

    """
    
    
    if type(fnames) == str:
        fnames = [fnames]
    
    for day_end in days_diff:
        for day_start in days_diff:
            if day_end < day_start:
                continue
            
            for fname in fnames:
        
                ds = xr.open_dataset(fpath + fname).load()
                
                df = add_var_from_nc(ds = ds, 
                                     variable = variable, 
                                     crs_ds = crs_4326, 
                                     df = df, 
                                     column = "{}_{}_{}".format(
                                                        column,
                                                        day_start,
                                                        day_end), 
                                     crs_df = crs_proj,
                                     days_off_start = day_start,
                                     days_off_end = day_end,
                                     function = function,
                                     offset = offset,
                                     init_dist = init_dist, 
                                     max_dist = max_dist,
                                     increment = increment)
                
    return df

def add_var_from_nc(ds, variable, crs_ds, df, column, crs_df, 
                    days_off_start = 0, days_off_end = 0,
                    function = "mean", offset = 0, 
                    init_dist = 500, max_dist = 5000, increment = 0):
    """
    Function to aggregate environmental data using a moving window around each
    in-situ measurement.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the environmental data.
    variable : str
        Name of the variable in 'ds' which shall be extracted.
    crs_ds : CRS
        CRS of 'ds'.
    df : pd.DataFrame
        DataFrame containg location and datetime of all samples for which the
        environmental data shall be extracted. Expects columns names "date",
        "east" and "north" (with UTM coordinates) or "date", "lat", "lon" (if 
        crs_proj is WGS84).
    column : str
        Name of the column which will be added to df and populated with the 
        data.
    crs_df : CRS
        CRS of the data in df..
    days_off_start : int, optional
        Window start (in days before measurement). The default is 0.
    days_off_end : int, optional
        Window end (in days before measurement). The default is 0.
    function : str, optional
        Function which will be used to aggregate data in the moving window. 
        Options are ["mean", "median", "max", "min", "sum" and "linreg"]. If 
        "linreg" is chosen, the slope of the linear regression is returned. The
        default is "mean".
    offset : int, optional
        Offset which will be added to data values (useful to convert K to °C). 
        The default is 0.
    init_dist : TYPE, optional
        Initial search radius around the Vibrio sampling location in the units
        of 'crs_df' (i.e. meters if in a UTM projection). The default is 500.
    max_dist : TYPE, optional
        Maximum search radius around the Vibrio sampling location in the units
        of 'crs_df' (i.e. meters if in a UTM projection). The default is 5000.
    increment : int, optional
        Increment by whih the search radius is increased in each iteration in 
        the units of 'crs_df' (i.e. meters if in a UTM projection). The default
        is 0.

    Returns
    -------
    pd.DataFrame
        Dataframe containing both Vibrio masurements and aggregated 
        environmental data for the requested seach window.
    """
    
    print("Adding data into column {}.".format(column))
    
    if "latitude" in ds:
        ds = ds.rename({"latitude": "lat", 
                        "longitude": "lon"})
        
    ds = ds.sortby(ds.time)        
    
    def agg_array(array, function = "mean"):
        if function == "mean":
            value = np.nanmean(array)
        elif function == "median":
            value = np.nanmedian(array)
        elif function == "max":
            value = np.nanmax(array)
        elif function == "min":
            value = np.nanmin(array)
        elif function == "sum":
            value = np.nansum(array)
        elif function == "linreg":
            # continue only if there are at least two measurements of which not
            # all are NaNs
            if array.size > 1 and not np.isnan(array).all():
                reg = LinearRegression()
                reg.fit(np.arange(len(array)).reshape(-1, 1), array)
                value = reg.coef_
            else:
                value = np.nan
        else:
            raise TypeError("Function {} not supported".format(function))
        return value                   

    if column not in df.columns:
        df[column] = np.nan
          
    # ignore numpy warnings for mean of empty slices
    warnings.filterwarnings(action='ignore', 
                            message='Mean of empty slice')
    warnings.filterwarnings(action='ignore', 
                            message='All-NaN slice encountered')
    
    for i, row in tqdm(df.iterrows(), total=df.shape[0]):
        try:
            ds_time = ds.sel(time = slice(
                        str(row.date - pd.DateOffset(days_off_end))[:10],
                        str(row.date - pd.DateOffset(days_off_start))[:10]))
            
            # check if there is data in the selection of the dataset
            if ("time" not in ds_time.sizes) or \
               ("time" in ds_time.sizes and ds_time.sizes["time"] > 0):
                
                if crs_ds == CRS("WGS84"):
                    value = agg_array(array = ds_time.sel(
                                        lat = row.lat,
                                        lon = row.lon,
                                        method = "nearest",
                                                  )[variable].mean().values,
                                     function = function)
                
                elif crs_ds == CRS("+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 "
                                   "+lat_2=63 +no_defs +R=6.371e+06"):
                    
                    transformer = Transformer.from_crs(crs_df, crs_ds)
                    
                    x, y = transformer.transform(row.east,
                                                 row.north)
                    
                    value = agg_array(array = ds_time.sel(
                                           x = x,
                                           y = y,
                                           method = "nearest",
                                           )[variable].values,
                                         function = function)
                
                else:
                    raise ValueError("Dataset has unexpected CRS!")
                
                # search for valid data in the area around the sampling 
                # location, if required
                if init_dist > 0:
                    
                    dist = init_dist
                    
                    while np.isnan(value):
                        
                        transformer = Transformer.from_crs(crs_df, crs_ds)
                        lat_min, lon_min = transformer.transform(row.east-dist, 
                                                                row.north-dist)
                        
                        lat_max, lon_max = transformer.transform(row.east+dist, 
                                                                row.north+dist)
                        
                        value = agg_array(array = ds_time.sel(
                                               lat = slice(lat_min, lat_max),
                                               lon = slice(lon_min, lon_max)
                                               )[variable].mean(["lat", "lon"]
                                                                ).values,
                                          function = function)
                
                        dist += increment
                        
                        if dist > max_dist:
                            break
                        
                if not np.isnan(value):
                    df.loc[i, column] = value + offset
        
        except KeyError:
            continue
    
    return df

if __name__ == "__main__":

    raw_df = pd.read_csv("N:/data/vibrio/cleaned/"
                            "vibrio_vuln_cleaned_sh_mv_categorical.csv", 
                            sep = "\t",
                            parse_dates = ["date"])


    #%% add mean SST from CMEMS data
    
    vib_mean_sst = raw_df.copy()
    
    fpath = "N:/data/SST_CMEMS/"
    fnames = ["SST_DMI_BAL_SST_L4_REP_OBSERVATIONS_010_016_2008-2019.nc",
              "SST_DMI-BALTIC-SST-L4-NRT-OBS_FULL_TIME_SERIE_2016-2021.nc"]
    
    vib_mean_sst = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = "analysed_sst", 
                                        df = raw_df.copy(), 
                                        column = "sst_cmems_mean", 
                                        days_diff = np.arange(0, 31), 
                                        function = "mean",
                                        crs_4326 = CRS("WGS84"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = -273.15,
                                        init_dist = 500, 
                                        max_dist = 5000, 
                                        increment = 500)
    
    vib_mean_sst.to_csv("N:/data/merged_vibrio_env/vibrio_mean_sst_cmems.csv", 
                        sep = "\t")                    
    
    #%% add max SST from CMEMS data
    
    fpath = "N:/data/SST_CMEMS/"
    fnames = ["SST_DMI_BAL_SST_L4_REP_OBSERVATIONS_010_016_2008-2019.nc",
              "SST_DMI-BALTIC-SST-L4-NRT-OBS_FULL_TIME_SERIE_2016-2021.nc"]
    
    vib_max_sst = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = "analysed_sst", 
                                        df = raw_df.copy(), 
                                        column = "sst_cmems_max", 
                                        days_diff = np.arange(30, 310, 10), 
                                        function = "max",
                                        crs_4326 = CRS("WGS84"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = -273.15,
                                        init_dist = 500, 
                                        max_dist = 5000, 
                                        increment = 500)
    
    vib_max_sst.to_csv("N:/data/merged_vibrio_env/vibrio_max_sst_cmems.csv", 
                        sep = "\t")
    
    #%% add SST Linear Regression from CMEMS data
    
    fpath = "N:/data/SST_CMEMS/"
    fnames = ["SST_DMI_BAL_SST_L4_REP_OBSERVATIONS_010_016_2008-2019.nc",
              "SST_DMI-BALTIC-SST-L4-NRT-OBS_FULL_TIME_SERIE_2016-2021.nc"]
    
    vib_linreg_sst = calc_env_data_offset(fpath = fpath, 
                                          fnames = fnames, 
                                          variable = "analysed_sst", 
                                          df = raw_df.copy(), 
                                          column = "sst_cmems_linreg", 
                                          days_diff = np.arange(0,31), 
                                          function = "linreg",
                                          crs_4326 = CRS("WGS84"), 
                                          crs_proj = CRS("EPSG:4647"), 
                                          offset = -273.15,
                                          init_dist = 500, 
                                          max_dist = 5000, 
                                          increment = 500)
    
    vib_linreg_sst.to_csv("N:/data/merged_vibrio_env/vibrio_sst_linreg_cmems.csv", 
                        sep = "\t")  
    
    
    
    #%% add mean SST from GETM data
    
    fpath = "D:/Vibrio/data/GETM/"
    fnames = ["WB200m.suf.2008-2021.nc"]
    
    vib_mean_sst = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = "temp", 
                                        df = raw_df.copy(), 
                                        column = "sst_mean", 
                                        days_diff = np.arange(0, 31), 
                                        function = "mean",
                                        crs_4326 = CRS("WGS84"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = 0,
                                        init_dist = 100, 
                                        max_dist = 300, 
                                        increment = 100)
    
    vib_mean_sst.to_csv("N:/data/merged_vibrio_env/vibrio_mean_sst_test.csv", 
                        sep = "\t")                    
    
    #%% add max SST from GETM data
    
    fpath = "D:/Vibrio/data/GETM/"
    fnames = ["WB200m.suf.2008-2021.nc"]
    
    vib_max_sst = calc_env_data_offset(fpath = fpath, 
                                       fnames = fnames, 
                                       variable = "temp", 
                                       df = raw_df.copy(), 
                                       column = "sst_max", 
                                       days_diff = np.arange(30, 310, 10), 
                                       function = "max",
                                       crs_4326 = CRS("WGS84"), 
                                       crs_proj = CRS("EPSG:4647"), 
                                       offset = 0,
                                       init_dist = 100, 
                                       max_dist = 300, 
                                       increment = 100)
    
    vib_max_sst.to_csv("N:/data/merged_vibrio_env/vibrio_max_sst_test.csv", 
                        sep = "\t")
    
    #%% add SST Linear Regression from GETM data
    
    fpath = "D:/Vibrio/data/GETM/"
    fnames = ["WB200m.suf.2008-2021.nc"]
    
    vib_linreg_sst = calc_env_data_offset(fpath = fpath, 
                                          fnames = fnames, 
                                          variable = "temp", 
                                          df = raw_df.copy(), 
                                          column = "sst_linreg", 
                                          days_diff = np.arange(0,31), 
                                          function = "linreg",
                                          crs_4326 = CRS("WGS84"), 
                                          crs_proj = CRS("EPSG:4647"), 
                                          offset = 0,
                                          init_dist = 100, 
                                          max_dist = 300, 
                                          increment = 100)
    
    vib_linreg_sst.to_csv("N:/data/merged_vibrio_env/vibrio_sst_linreg_test.csv", 
                        sep = "\t")  
    
    #%% add mean temperature from MetPP
    
    fpath = "N:/data/meteo_data/"
    fnames = ["MetPPArchive2_daily_avg_2013-2021.nc"]
    
    vib_mean_at = calc_env_data_offset(fpath = fpath, 
                                       fnames = fnames, 
                                       variable = "air_temperature_2m_mean", 
                                       df = raw_df.copy(), 
                                       column = "AT_mean", 
                                       days_diff = np.arange(0,31), 
                                       function = "mean",
                                       crs_4326 = CRS("+proj=lcc +lat_0=63 "
                                                      "+lon_0=15 +lat_1=63 "
                                                      "+lat_2=63 +no_defs "
                                                      "+R=6.371e+06"), 
                                       crs_proj = CRS("EPSG:4647"), 
                                       offset = -273.15,
                                       init_dist = 0, 
                                       max_dist = 0, 
                                       increment = 0)
    
    vib_mean_at.to_csv("N:/data/merged_vibrio_env/vibrio_mean_at.csv", 
                        sep = "\t")
    
    #%% add mean delta temperature from MetPP
    
    fpath = "N:/data/meteo_data/"
    fnames = ["MetPPArchive2_daily_avg_2013-2021.nc"]
    
    vib_mean_at = calc_env_data_offset(fpath = fpath, 
                                       fnames = fnames, 
                                       variable = "air_temperature_2m_mean", 
                                       df = raw_df.copy(), 
                                       column = "AT_linreg", 
                                       days_diff = np.arange(0,31), 
                                       function = "linreg",
                                       crs_4326 = CRS("+proj=lcc +lat_0=63 "
                                                      "+lon_0=15 +lat_1=63 "
                                                      "+lat_2=63 +no_defs "
                                                      "+R=6.371e+06"), 
                                       crs_proj = CRS("EPSG:4647"), 
                                       offset = -273.15,
                                       init_dist = 0, 
                                       max_dist = 0, 
                                       increment = 0)
    
    vib_mean_at.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_at.csv", 
                        sep = "\t")  
    
    
    #%% add irradiation from MetPP
    
    fpath = "N:/data/meteo_data/"
    fnames = ["MetPPArchive2_daily_avg_2013-2021.nc"]
    variable = "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time_sum"
    
    vib_mean_irr = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = variable, 
                                        df = raw_df.copy(), 
                                        column = "irrad_mean", 
                                        days_diff = np.arange(0,31), 
                                        function = "mean",
                                        crs_4326 = CRS("+proj=lcc +lat_0=63 "
                                                       "+lon_0=15 +lat_1=63 "
                                                       "+lat_2=63 +no_defs "
                                                       "+R=6.371e+06"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = 0,
                                        init_dist = 0, 
                                        max_dist = 0, 
                                        increment = 0)
    
    vib_mean_irr.to_csv("N:/data/merged_vibrio_env/vibrio_mean_irr.csv", 
                        sep = "\t")
    
    #%% add delta of irradiation from MetPP
    
    fpath = "N:/data/meteo_data/"
    fnames = ["MetPPArchive2_daily_avg_2013-2021.nc"]
    variable = "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time_sum"
    
    vib_mean_irr = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = variable, 
                                        df = raw_df.copy(), 
                                        column = "irrad_linreg", 
                                        days_diff = np.arange(0,31), 
                                        function = "linreg",
                                        crs_4326 = CRS("+proj=lcc +lat_0=63 "
                                                       "+lon_0=15 +lat_1=63 "
                                                       "+lat_2=63 +no_defs "
                                                       "+R=6.371e+06"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = 0,
                                        init_dist = 0, 
                                        max_dist = 0, 
                                        increment = 0)
    
    vib_mean_irr.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_irr.csv", 
                        sep = "\t") 
    
    #%% add wind speed from MetPP
    
    fpath = "N:/data/meteo_data/"
    fnames = ["MetPPArchive2_daily_avg_2013-2021.nc"]
    
    vib_mean_wind = calc_env_data_offset(fpath = fpath, 
                                         fnames = fnames, 
                                         variable = "wind_speed_10m_mean", 
                                         df = raw_df.copy(), 
                                         column = "ws_mean", 
                                         days_diff = np.arange(0,31), 
                                         function = "mean",
                                         crs_4326 = CRS("+proj=lcc +lat_0=63 "
                                                        "+lon_0=15 +lat_1=63 "
                                                        "+lat_2=63 +no_defs "
                                                        "+R=6.371e+06"), 
                                         crs_proj = CRS("EPSG:4647"), 
                                         offset = 0,
                                         init_dist = 0, 
                                         max_dist = 0, 
                                         increment = 0)
    
    vib_mean_wind.to_csv("N:/data/merged_vibrio_env/vibrio_mean_wind_speed.csv", 
                        sep = "\t")
    
    #%% add delta of wind speed from MetPP
    
    fpath = "N:/data/meteo_data/"
    fnames = ["MetPPArchive2_daily_avg_2013-2021.nc"]
    
    vib_linreg_wind = calc_env_data_offset(fpath = fpath, 
                                           fnames = fnames, 
                                           variable = "wind_speed_10m_mean", 
                                           df = raw_df.copy(), 
                                           column = "ws_linreg", 
                                           days_diff = np.arange(0,31), 
                                           function = "linreg",
                                           crs_4326 = CRS("+proj=lcc +lat_0=63 "
                                                          "+lon_0=15 +lat_1=63 "
                                                          "+lat_2=63 +no_defs "
                                                          "+R=6.371e+06"), 
                                           crs_proj = CRS("EPSG:4647"), 
                                           offset = 0,
                                           init_dist = 0, 
                                           max_dist = 0, 
                                           increment = 0)
    
    vib_linreg_wind.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_wind_speed.csv", 
                            sep = "\t") 
    
    #%% add precipitation from MetPP
    
    fpath = "N:/data/meteo_data/"
    fnames = ["MetPPArchive2_daily_avg_2013-2021.nc"]
    
    vib_mean_prec = calc_env_data_offset(fpath = fpath, 
                                         fnames = fnames, 
                                         variable = "precipitation_amount_sum", 
                                         df = raw_df.copy(), 
                                         column = "prec_mean", 
                                         days_diff = np.arange(0,31), 
                                         function = "mean",
                                         crs_4326 = CRS("+proj=lcc +lat_0=63 "
                                                        "+lon_0=15 +lat_1=63 "
                                                        "+lat_2=63 +no_defs "
                                                        "+R=6.371e+06"), 
                                         crs_proj = CRS("EPSG:4647"), 
                                         offset = 0,
                                         init_dist = 0, 
                                         max_dist = 0, 
                                         increment = 0)
    
    vib_mean_prec.to_csv("N:/data/merged_vibrio_env/vibrio_mean_precipitation.csv", 
                        sep = "\t")
    
    #%% add delta of wind speed from MetPP
    
    fpath = "N:/data/meteo_data/"
    fnames = ["MetPPArchive2_daily_avg_2013-2021.nc"]
    
    vib_linreg_prec = calc_env_data_offset(fpath = fpath, 
                                           fnames = fnames, 
                                           variable = "precipitation_amount_sum", 
                                           df = raw_df.copy(), 
                                           column = "prec_linreg", 
                                           days_diff = np.arange(0,31), 
                                           function = "linreg",
                                           crs_4326 = CRS("+proj=lcc +lat_0=63 "
                                                          "+lon_0=15 +lat_1=63 "
                                                          "+lat_2=63 +no_defs "
                                                          "+R=6.371e+06"), 
                                           crs_proj = CRS("EPSG:4647"), 
                                           offset = 0,
                                           init_dist = 0, 
                                           max_dist = 0, 
                                           increment = 0)
    
    vib_linreg_prec.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_precipitation.csv", 
                            sep = "\t") 
    
    
    #%% add mean SSS from GETM data
    
    fpath = "N:/data/GETM/"
    fnames = ["WB200m.suf.2008-2021.nc"]
    
    vib_mean_sal = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = "salt", 
                                        df = raw_df.copy(), 
                                        column = "sal_mean", 
                                        days_diff = np.arange(0,31), 
                                        function = "mean",
                                        crs_4326 = CRS("WGS84"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = 0,
                                        init_dist = 100, 
                                        max_dist = 300, 
                                        increment = 100)
    
    vib_mean_sal.to_csv("N:/data/merged_vibrio_env/vibrio_mean_sal.csv", 
                        sep = "\t")
    
    #%% add delta SSS from GETM data
    
    fpath = "D:/Vibrio/data/GETM/"
    fnames = ["WB200m.suf.2008-2021 - Copy.nc"]
    
    vib_delta_sal = calc_env_data_offset(fpath = fpath, 
                                         fnames = fnames, 
                                         variable = "salt", 
                                         df = raw_df.copy(), 
                                         column = "sal_linreg", 
                                         days_diff = np.arange(0,31), 
                                         function = "linreg",
                                         crs_4326 = CRS("WGS84"), 
                                         crs_proj = CRS("EPSG:4647"), 
                                         offset = 0,
                                         init_dist = 100, 
                                         max_dist = 300, 
                                         increment = 100)
    
    vib_delta_sal.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_sal_test.csv", 
                          sep = "\t")
        
    #%% add Chl from CMEMS
    
    fpath = "N:/data/nutrients/"
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_mean_chl = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = "chl", 
                                        df = raw_df.copy(), 
                                        column = "chl_mean", 
                                        days_diff = np.arange(0,31), 
                                        function = "mean",
                                        crs_4326 = CRS("WGS84"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = 0,
                                        init_dist = 2000, 
                                        max_dist = 6000, 
                                        increment = 1000)
    
    vib_mean_chl.to_csv("N:/data/merged_vibrio_env/vibrio_mean_chl.csv", 
                        sep = "\t")
    
    #%% add delta Chl from CMEMS
    
    fpath = "N:/data/nutrients/"
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_linreg_chl = calc_env_data_offset(fpath = fpath,
                                          fnames = fnames,
                                          variable = "chl",
                                          df = raw_df.copy(),
                                          column = "chl_linreg",
                                          days_diff = np.arange(0,31),
                                          function = "linreg",
                                          crs_4326 = CRS("WGS84"),
                                          crs_proj = CRS("EPSG:4647"),
                                          offset = 0,
                                          init_dist = 2000, 
                                          max_dist = 6000, 
                                          increment = 1000)
    
    vib_linreg_chl.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_chl.csv", 
                        sep = "\t") 
    
    #%% add O2 from CMEMS
    
    fpath = "N:/data/nutrients/"
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_mean_o2 = calc_env_data_offset(fpath = fpath, 
                                       fnames = fnames, 
                                       variable = "o2", 
                                       df = raw_df.copy(), 
                                       column = "o2_mean", 
                                       days_diff = np.arange(0,31), 
                                       function = "mean",
                                       crs_4326 = CRS("WGS84"), 
                                       crs_proj = CRS("EPSG:4647"), 
                                       offset = 0,
                                       init_dist = 2000, 
                                       max_dist = 6000, 
                                       increment = 1000)
    
    vib_mean_o2.to_csv("N:/data/merged_vibrio_env/vibrio_mean_o2.csv", 
                        sep = "\t") 
    
    #%% add delta O2 from CMEMS
    
    fpath = "N:/data/nutrients/"
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_linreg_o2 = calc_env_data_offset(fpath = fpath, 
                                         fnames = fnames, 
                                         variable = "o2", 
                                         df = raw_df.copy(), 
                                         column = "o2_linreg", 
                                         days_diff = np.arange(0,31), 
                                         function = "linreg",
                                         crs_4326 = CRS("WGS84"), 
                                         crs_proj = CRS("EPSG:4647"), 
                                         offset = 0,
                                         init_dist = 2000, 
                                         max_dist = 6000, 
                                         increment = 1000)
    
    vib_linreg_o2.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_o2.csv", 
                        sep = "\t") 
    
    #%% add NH4 from CMEMS
    
    fpath = "N:/data/nutrients/"
    
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_mean_nh4 = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = "nh4", 
                                        df = raw_df.copy(), 
                                        column = "nh4_mean", 
                                        days_diff = np.arange(0,31), 
                                        function = "mean",
                                        crs_4326 = CRS("WGS84"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = 0,
                                        init_dist = 2000, 
                                        max_dist = 6000, 
                                        increment = 1000)
    
    vib_mean_nh4.to_csv("N:/data/merged_vibrio_env/vibrio_mean_nh4.csv", 
                        sep = "\t") 
    
    #%% add delta NH4 from CMEMS
    
    fpath = "N:/data/nutrients/"
    
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_linreg_nh4 = calc_env_data_offset(fpath = fpath, 
                                          fnames = fnames, 
                                          variable = "nh4", 
                                          df = raw_df.copy(), 
                                          column = "nh4_linreg", 
                                          days_diff = np.arange(0,31), 
                                          function = "linreg",
                                          crs_4326 = CRS("WGS84"), 
                                          crs_proj = CRS("EPSG:4647"), 
                                          offset = 0,
                                          init_dist = 2000, 
                                          max_dist = 6000, 
                                          increment = 1000)
    
    vib_linreg_nh4.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_nh4.csv", 
                        sep = "\t")
    
    #%% add NO3 from CMEMS
    
    fpath = "N:/data/nutrients/"
    
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_mean_no3 = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = "no3", 
                                        df = raw_df.copy(), 
                                        column = "no3_mean", 
                                        days_diff = np.arange(0,31), 
                                        function = "mean",
                                        crs_4326 = CRS("WGS84"), 
                                        crs_proj = CRS("EPSG:4647"), 
                                        offset = 0,
                                        init_dist = 2000, 
                                        max_dist = 6000, 
                                        increment = 1000)
    
    vib_mean_no3.to_csv("N:/data/merged_vibrio_env/vibrio_mean_no3.csv", 
                        sep = "\t") 
    
    #%% add delta NO3 from CMEMS
    
    fpath = "N:/data/nutrients/"
    
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_linreg_no3 = calc_env_data_offset(fpath = fpath, 
                                          fnames = fnames, 
                                          variable = "no3", 
                                          df = raw_df.copy(), 
                                          column = "no3_linreg", 
                                          days_diff = np.arange(0,31), 
                                          function = "linreg",
                                          crs_4326 = CRS("WGS84"), 
                                          crs_proj = CRS("EPSG:4647"), 
                                          offset = 0,
                                          init_dist = 2000, 
                                          max_dist = 6000, 
                                          increment = 1000)
    
    vib_linreg_no3.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_no3.csv", 
                        sep = "\t") 
    
    #%% add po4 from CMEMS
    
    fpath = "N:/data/nutrients/"
    
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_mean_po4 = calc_env_data_offset(fpath = fpath, 
                                        fnames = fnames, 
                                        variable = "po4", 
                                        df = raw_df.copy(), 
                                        column = "po4_mean", 
                                        days_diff = np.arange(0,31), 
                                        function = "mean",
                                        crs_4326 = CRS("WGS84"),
                                        crs_proj = CRS("EPSG:4647"),
                                        offset = 0,
                                        init_dist = 2000, 
                                        max_dist = 6000, 
                                        increment = 1000)
    
    vib_mean_po4.to_csv("N:/data/merged_vibrio_env/vibrio_mean_po4.csv", 
                        sep = "\t") 
    
    #%% add delta po4 from CMEMS
    
    fpath = "N:/data/nutrients/"
    
    fnames = ["CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc"]
    
    vib_linreg_po4 = calc_env_data_offset(fpath = fpath, 
                                          fnames = fnames, 
                                          variable = "po4", 
                                          df = raw_df.copy(), 
                                          column = "po4_linreg", 
                                          days_diff = np.arange(0,31), 
                                          function = "linreg",
                                          crs_4326 = CRS("WGS84"), 
                                          crs_proj = CRS("EPSG:4647"), 
                                          offset = 0,
                                          init_dist = 2000, 
                                          max_dist = 6000, 
                                          increment = 1000)
    
    vib_linreg_po4.to_csv("N:/data/merged_vibrio_env/vibrio_linreg_po4.csv", 
                        sep = "\t") 


