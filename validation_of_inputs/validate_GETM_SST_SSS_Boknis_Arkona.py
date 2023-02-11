# -*- coding: utf-8 -*-
"""
Carries out the validation for Boknis Eck, Kühlungsborn and Arkona Basin 
stations using the match_up class class from `model_match_up`. Finally, the 
results are plotted in a sactter plot (Figure 3).
"""

import pandas as pd
import geopandas as gpd
from dask.diagnostics import ProgressBar
import xarray as xr
import matplotlib.pyplot as plt
from model_match_up import match_up

#%% This takes long. Run only when necessary!
# open and prepare model data
with ProgressBar():
        xds = xr.open_mfdataset('N:/data/GETM/*.nc4', chunks = "auto")
        xds = xds.rename({"latc": "lat", "lonc": "lon"})
        #xds = xds.sel( time=slice('2018-01-01', '2019-01-01')).load()
        lats = xr.DataArray([54.156053, 54.5295, 54.6967, 54.8833 ],
                            dims='lat')
        lons = xr.DataArray([10.0393, 11.763174, 12.705, 13.8667], 
                            dims='lon')
        xds = xds.sel(lon=lons, lat=lats, method='nearest').load()

xds.rio.write_crs(4326, inplace=True)

#%% define some parameters for the plotting 
ms_bok = 1 # markersize for Boknis Eck
al_bok = 1 # alpha for Boknis Eck
ms = 0.1 # markersize for Arkona
al = 0.5 # alpha for Arkona
c_bok = "tab:blue"#66c2a5"#tab:blue"
c_kul = "tab:gray"#fc8d62"#tab:pink"
c_ark = "tab:orange"#8da0cb"#tab:gray"

#%% open and prepare Bokniseck in-situ data for SSS
in_situ_fpath = "N:/data/in-situ/boknis-eck/boknis_eck_data_export_2023-01-13_sal.csv"
data_col = "ctdsal"
time_col = "Time"

df_bok = pd.read_csv(in_situ_fpath, parse_dates = [time_col])
df_bok = df_bok[df_bok["ctdsal_flag"] == 1]
df_bok = (df_bok
          .drop(["Name"], axis = 1)
          .groupby("Label")
          .mean(numeric_only=False))
df_bok = df_bok.rename(columns={"ctdsal": "sal"})
df_bok["station"] = c_bok
df_bok["alpha"] = al_bok
df_bok["markersize"] = ms_bok 


#%% open and prepare Arkona in-situ data for SSS
in_situ_fpath = "N:/data/in-situ/arkona_basin/schuett.20230119120856495.arkona_basin_psal_2m.nc"
data_col = "sal"
time_col = "Time"

xds_in = xr.open_dataset(in_situ_fpath)
df_in = xds_in.to_dataframe()
qc_df = xds_in["PSAL_QC.MicroCAT"].astype(int).to_dataframe()
df_in["PSAL_QC.MicroCAT"] = qc_df["PSAL_QC.MicroCAT"]
df_in = df_in[(df_in["PSAL_QC.MicroCAT"] == 1) | 
              (df_in["PSAL_QC.MicroCAT"] == 2)]
xds_in = df_in.to_xarray()
xds_in = xds_in.resample(TIME='1D').mean('TIME')
df = xds_in.to_dataframe().reset_index()
df["TIME"] = pd.to_datetime(df["TIME"])
df = df.rename(columns={"LONGITUDE": "Longitude", 
                        "LATITUDE": "Latitude", 
                        "PSAL.MicroCAT": "sal", 
                        "TIME": "Time"})
df["station"] = c_ark
df["alpha"] = al
df["markersize"] = ms

# #%% open and prepare Darss Sill in-situ data for SSS
# df_ds = pd.read_excel("N:/data/in-situ/darss_sill/odin2_2023-02-08_095537_cleaned.xlsx")
# df_ds = (df_ds.set_index("Time", drop = True)
#          .resample('D')
#          .mean(numeric_only = True)
#          .reset_index()
#          .dropna())
# df_ds = df_ds[["Time", "PSAL5DSD", "TEMP5STD"]].rename({"TEMP5STD": "temp",
#                                                         "PSAL5DSD": "sal"}, 
#                                                        axis = 1)
# df_ds["Longitude"] = 12.705
# df_ds["Latitude"] = 54.6967
# df_ds["station"] = "tab:red"
# df_ds["alpha"] = al
# df_ds["markersize"] = ms
 
#%%  run match-up pipeline for SSS


# concatenate both dataframes
jdf = pd.concat([df,df_bok, #df_ds
                 ]).reset_index(drop = True)
jgdf = gpd.GeoDataFrame(jdf, 
                        geometry = gpd.points_from_xy(jdf["Longitude"], 
                                                      jdf["Latitude"]),
                        crs = "4326").dropna(subset = ['Latitude',"Longitude"])

xds_var = "salt" 
# generate match-up validation
sss = match_up(gdf = jgdf,
                time_col = time_col, 
                var_dict = {data_col: xds_var}, 
                xds = xds,
                station = "station",
                alpha = "alpha",
                markersize = "markersize",
                max_dist = 200,
                verbose = 1)

sssdf = sss.data["salt"]
sssdf.to_csv("N:/data/in-situ/sss_matchup.csv")


#%% do the same of SST

# open and prepare Arkona in-situ data for SST
in_situ_fpath = "N:/data/in-situ/arkona_basin/schuett.20230119163626580.arkona_basin_temp_05m.nc"

xds_in = xr.open_dataset(in_situ_fpath)
df_in = xds_in.to_dataframe()
qc_df = xds_in["TEMP_QC.DWR"].astype(int).to_dataframe()
df_in["TEMP_QC.DWR"] = qc_df["TEMP_QC.DWR"]
df_in = df_in[(df_in["TEMP_QC.DWR"] == 1) | (df_in["TEMP_QC.DWR"] == 2)]
xds_in = df_in.to_xarray()
xds_in = xds_in.resample(TIME='1D').mean('TIME')
df = xds_in.to_dataframe().reset_index()
df["TIME"] = pd.to_datetime(df["TIME"])
df = df.rename(columns={"LONGITUDE": "Longitude", 
                        "LATITUDE": "Latitude", 
                        "TEMP.DWR": "temp", 
                        "TIME": "Time"})
df["station"] = c_ark
df["alpha"] = al
df["markersize"] = ms

xds_var = "temp"

#%% open and prepare SST data from Kühlungsborn
in_situ_fpath = "N:/data/in-situ/kuehlungsborn/sst.kuehlungsborn.nc"
xds_in = xr.open_dataset(in_situ_fpath)
df_kul = xds_in.to_dataframe()
df_kul = df_kul[df_kul["TEMP_QC"] <= 2]
df_kul = (df_kul.set_index("time", drop = True)
         .resample('D')
         .mean(numeric_only = True)
         .reset_index()
         .dropna())
df_kul = df_kul[["TEMP", "time"]].rename({"TEMP": "temp",
                                          "time": "Time"}, axis = 1)
df_kul["Longitude"] = 11.763174
df_kul["Latitude"] = 54.156053
df_kul["station"] = c_kul
df_kul["alpha"] = al_bok
df_kul["markersize"] = ms


#%% open and prepare Bokniseck in-situ data for SST
in_situ_fpath = "N:/data/in-situ/boknis-eck/boknis_eck_data_export_2023-01-13_temp.csv"
data_col = "temp"
time_col = "Time"

df_bok = pd.read_csv(in_situ_fpath, parse_dates = [time_col])
df_bok = df_bok[df_bok["ctdtmp_flag"] == 1]
df_bok = df_bok.drop(["Name"], axis = 1).groupby("Label").mean(numeric_only=False)
df_bok = df_bok.rename(columns={"ctdtmp": "temp"})
df_bok["station"] = c_bok
df_bok["alpha"] = al_bok
df_bok["markersize"] = ms_bok

#%% join dataframes from Boknis Eck and Arkona Basin
jdf = pd.concat([df_kul,
                 df, #df_ds, 
                 df_bok])
jgdf = gpd.GeoDataFrame(jdf, 
                        geometry = gpd.points_from_xy(jdf["Longitude"], 
                                                      jdf["Latitude"]),
                        crs = "4326").dropna(subset = ['Latitude',"Longitude"])


#%% run the match-up validation for SST
sst = match_up(gdf = jgdf,
               time_col = time_col, 
               var_dict = {data_col: xds_var}, 
               xds = xds,
               station = "station",
               alpha = "alpha",
               markersize = "markersize",
               # xds_var = xds_var,
               max_dist = 200,
               verbose = 1)

sstdf = sst.data["temp"]
sstdf.to_csv("N:/data/in-situ/sst_matchup.csv")


#%% plot results
# import seaborn as sns

fig, [ax1, ax2] = plt.subplots(1,2, figsize = (8,4))

sst.plot(ax = ax1,
         key = "temp",
         plot_stats = True)

# sns.kdeplot(ax = ax1,
#             data = sstdf, x="observed", y="modelled",
#             levels = [0.25, 0.5, 0.75],
#             color = "k",
#             linewidths = 0.5)

ax1.set_title("SST [°C]")

sss.plot(ax = ax2,
         key = "salt",
         ylabel = "",
         plot_stats = True)

# sns.kdeplot(ax = ax2,
#             data = sssdf, x="observed", y="modelled",
#             levels = [0.25, 0.5, 0.75],
#             color = "k",
#             linewidths = 0.5)

ax2.set_title("SSS [PSU]")


plt.savefig("N:/plots/GETM-Boknis_Eck_Arkona_SST_SSS.jpeg",
            bbox_inches = "tight",
            dpi = 300)


