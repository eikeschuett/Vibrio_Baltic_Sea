# -*- coding: utf-8 -*-
"""
Compiles the data of all environmental parameters and exports them. 
Additionally, a csv containing only the optimal time lags will be exported.
"""

import os
import pandas as pd

def prep_df(df):
    df = df.drop(['Unnamed: 0.1', 'Unnamed: 0'], axis = 1)
    df["date"] = [pd.to_datetime(x[0:10]) for x in df["date"]]
    df = df[df["date"] < "2019-09-01"]
    df = df.drop_duplicates()
    df = df.reset_index(drop = True)
    return df

#import data and prepare for analysis
os.chdir("N:/data/merged_vibrio_env/")

files = [x for x in os.listdir() if x.startswith(("vibrio_mean_",
                                                  "vibrio_linreg_",
                                                  "vibrio_max",
                                                  "vibrio_sst_linreg"))]

on_cols = ['date', 'station_id', 'station_name', 'species', 'quantity',
            'temperature', 'east', 'north', 'lat', 'lon', 'old_id', 'comments',
            'air temp', 'salinity']

for file in files:
    if file == files[0]:
        data = pd.read_csv(file, sep = "\t")
        print(len(data))
        data = prep_df(df = data)
        
    else:
        df = pd.read_csv(file, sep = "\t")
        print(len(df))
        df = prep_df(df = df)
        
        data = pd.merge(data, df, on = on_cols)
        
data = data.iloc[data.drop(["old_id"], axis = 1).drop_duplicates(
    ).index].reset_index(drop = True)

# final filter for stations with unreliable salinity
min_sal = data.groupby(["station_id"]).agg(min = ("sal_mean_0_0", "mean"))
sal_filter = min_sal[(min_sal["min"] < 2) | (min_sal["min"].isna())]
data = data[~data["station_id"].isin(sal_filter.index)]

# export full dataframe (i.e. all time lag combinations)
data.to_csv("vibrio_data_merged_sh_mv.csv", 
            index = False, 
            sep = "\t")

# select optimal time lags from the dataframe
data = data[["quantity", 
                    
             # floor
             "chl_mean_6_15", 
             "chl_linreg_3_5", 
             "o2_mean_0_0",
             "o2_linreg_5_16", 
             "no3_mean_7_25",
             "no3_linreg_0_24",
             "nh4_mean_28_28",
             "nh4_linreg_9_10",
             "po4_mean_22_22", 
             "po4_linreg_4_9", 
             "sst_mean_0_11", 
             "sst_linreg_0_30",
             "sst_max_180_180",
             "sal_mean_0_6", 
             "sal_linreg_11_26",
             "ASWDIFD_S_mean_27_28",
             "ASWDIFD_S_linreg_18_28",
             "T_2M_CL_mean_0_16",#
             "T_2M_CL_linreg_7_29",
             "tp_mean_6_18", 
             "tp_linreg_5_8",
             "ws_mean_mean_2_17",
             "ws_mean_linreg_4_8",             
             ]]

data = data.dropna()

data.to_csv("vibrio_data_merged_sh_mv_timelag.csv", 
            sep = "\t",
            index = False,
            )    