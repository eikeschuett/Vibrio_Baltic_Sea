
"""
Calculates Vibrio season length, apply Sen's slope and Men-Kendall trend 
analysis and finally plot the data (Figure 7).
"""

import numpy as np
import pandas as pd 
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
import sys
sys.path.append("your/path/to/the/repository/plots")
from Fig1_4_map_study_area import (create_map, add_aligned_colorbar, 
                                   remove_Schlei, add_subplot_char)
from cmocean import cm
from cmocean.tools import crop_by_percent
import os
import pymannkendall as mk
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix

def rolling_sst_sss(fpath, file):
    """
    Calculates a rolling mean of SST and SSS (with the time lags defined in 
    the correlation matrices, i.e. 0-6 days for SSS and 0-11 days for SST)
    """
    with xr.open_dataset(os.path.join(fpath, file), chunks = "auto") as xds:
        with ProgressBar():
            xds_2 = xds["salt"].rolling(time = 6, 
                                        min_periods = 6, 
                                        center = False).mean().to_dataset(
                                            name = "salt_0_6").compute()
            
        with ProgressBar():    
            xds_2["temp_0_11"] = xds["temp"].rolling(time = 11, 
                                                     min_periods = 11, 
                                                     center = False).mean(
                                                         ).compute()
        
        xds_2.to_netcdf(os.path.join(fpath, "rolling_mean/", file[:-4] + 
                                     "_mean.nc"))
        
        return file
    
def mannkendall(xi, yi):
    """
    Define a function to run the Mann-Kendall text using the PyMannKendall
    implementation (https://doi.org/10.21105/joss.01556). It includes the 
    calculation of Sen's slope.

    Parameters
    ----------
    xi : np.array or similar
        Values of the x-axis (in this case year).
    yi : np.array or similar
        Values of the y-axis (in this case SST).

    Returns
    -------
    trend : np.array
        significant trends and their direction.
    h : np.array
        boolean if H0 (no trend) or Ha (alternative hypothesis) are accepted.
    p : np.array
        significance of the trend.
    z : np.array
        Positive/negative z: data values tend to increase/decrease with time.
    Tau : np.array
        proportion of up-movements against time vs the proportion of down-
        movements.
    s : np.array
        Mann-Kendal's score.
    var_s : np.array
        Variance S.
    slope : np.array
        Theil-Sen estimator/slope.
    intercept : TYPE
        intercept of Kendall-Theil Robust Line (note: intercept is for year 
        1995, not for year 0).
    """

    trend, h, p, z, Tau, s, var_s, slope, intercept  = mk.original_test(yi)
    return trend, h, p, z, Tau, s, var_s, slope, intercept

def round_off_rating(number):
    """
    Round a number to the closest half integer.
    """
    return round(number * 2) / 2

def logistic_regression(df):

    df = df[["quantity", "sst_mean_0_11"]]
    
    df["positive"] = 0
    df["positive"][df["quantity"] > 0] = 1
    df["sst_round"] = round_off_rating(df['sst_mean_0_11'])
    
    df_grouped = df.groupby('sst_round').mean()
    
    model = LogisticRegression(solver='lbfgs', random_state=0)
    
    x = np.array(df["sst_mean_0_11"]).reshape(-1, 1)
    y = np.array(df["positive"])
    
    # fit the model
    model.fit(x, y)
    
    # get the score of the model (i.e. percent of correct classifications)
    # model.score(x, y)
    
    # print confusion matrix
    # print(confusion_matrix(y, model.predict(x)))
    # i.e. 230 true negatives, 80 false negatives, 88 false positives, 
    # 223 true positives
    
    # print(classification_report(y, model.predict(x)))
    
    return model, df_grouped

def plot_log_fig(model, df_grouped):
    #plot the fit for the appendix:
    # predict the probabilities for a range of SST to plot them
    x_proba = np.arange(7, 25, 0.1).reshape(-1,1)
    proba = [x[1] for x in model.predict_proba(x_proba)]
    
    # plot the logistic fit
    fig, ax = plt.subplots(figsize = (6,3))
    ax.plot(df_grouped.index, df_grouped.positive)
    
    # add vertical and horizontal lines at our SST threshold
    ax.axhline(y = 0.33, color = 'gray', linestyle = '-', linewidth = 0.5)
    ax.axvline(x = 17.6, color = 'gray', linestyle = '-', linewidth = 0.5)
    
    #plot the logistic fit
    ax.plot(x_proba, proba)
    
    ax.set_xlabel("$\overline{SST}_{0-11}$")
    ax.set_ylabel("Vibrio occurence probability [-]")
    
    ax.set_xlim([10, 24.9])
    ax.set_ylim([0, 1])
    
    plt.savefig("N:/plots/log_fit_vibrio_occurence.jpeg", 
                dpi = 300,
                bbox_inches = "tight")
    plt.show()


#%%

if __name__ == "__main__":
    
    # logistic fit to derive the SST threshold
    df = pd.read_csv(r"N:\data\merged_vibrio_env\vibrio_data_merged_sh_mv_timelag.csv",
                     sep = "\t")
    
    model, df_grouped = logistic_regression(df)
    
    # # probability of Vibrio occurence at 17.6 °C is ~ 0.33
    print(model.predict_proba(np.array([17.6]).reshape(-1,1)))
    
    plot_log_fig(model, df_grouped)
        
    # define the threshold for SST (17.6 °C is the 33 % occurence probability)
    thres = 17.6
    
    #%% calculate rolling means of SST and SSS
    
    fpath = "D:/Vibrio/data/GETM"
    files = [x for x in os.listdir(fpath) if x.endswith(".nc4")]
    
    print([rolling_sst_sss(fpath, file) for file in files])
    
    #%% calculate Vibrio season length based on the SST rolling mean
    
    xds = xr.open_mfdataset('D:/Vibrio/data/GETM/rolling_mean/*.nc', 
                            chunks = "auto")
    xds.rio.write_crs(ccrs.UTM(zone = 32), inplace=True)
    xds = xds.sel(time = slice("1995-01-01", "2021-12-31"))
    xds = remove_Schlei(xds,xdim = "lonc", ydim = "latc")
    
    # calculate Vibrio season length                    
    with ProgressBar():
        day_count = xds.temp_0_11 \
                        .where(xds.temp_0_11 >= thres) \
                        .groupby('time.year') \
                        .count(dim='time')\
                        .compute()
     
    #%% apply the Mann-Kendall trend test to our dataset
    # takes ~ 15 minutes on my PC
    
    trend, h, p, z, Tau, s, var_s, slope, intercept = xr.apply_ufunc(
        mannkendall,
        day_count['year'],
        day_count,
        input_core_dims=[['year'],['year']],
        output_core_dims=[[],[],[],[],[],[],[],[],[]],
        vectorize=True,
        dask='parallelized',
        )
    
    #%% create dataset from individual data arrays and save them
    xds_result = xr.Dataset({"day_count": day_count,
                             "slope": slope,
                             "p": p})
    
    xds_result.to_netcdf("N:/data/vibrio_season/season_trend_analysis_sst_17.6_w_lag.nc")
    
    #%% plot the results
    xds_result = xr.open_dataset("N:/data/vibrio_season/season_trend_analysis_sst_17.6_w_lag.nc").load()
    day_count = xds_result.day_count
    slope = xds_result.slope 
    p = xds_result.p
    
    extent = [525000, 850000, 5963000, 6091000]
    
    fig = plt.figure(figsize=(8,3.25))
    
    fig, ax = create_map(fig = fig,
                          sp_x = 2, 
                          sp_y = 2, 
                          sp_n = 1,
                          extent = extent,
                          bottom_labels = False,
                          gridticksize = 8)
    add_subplot_char(ax, "a", loc = "lower left")
    
    pc = day_count.sel(year = 1995).where(day_count.sel(year = 1995) > 0
                                              ).plot(ax = ax,
                                                add_colorbar=False,
                                                transform=ccrs.PlateCarree(),
                                                vmin = 45,
                                                vmax = 105,
                                                )
    ax.set_title("")
    
    cbar, ax_cb = add_aligned_colorbar(fig = fig, ax = ax, pc = pc)
    
    cbar.set_label(label = "length of ${V. vulnificus}$ \n "
            "season 1995 [days]",size = 8)
    
    
    fig, ax = create_map(fig = fig,
                          sp_x = 2, 
                          sp_y = 2, 
                          sp_n = 2,
                          extent = extent,
                          bottom_labels = False,
                          left_labels = False,
                          gridticksize = 8)
    add_subplot_char(ax, "b", loc = "lower left")
    
    pc = day_count.sel(year = 2021).where(day_count.sel(year = 2021)\
                                              ).plot(ax = ax,
                                                add_colorbar=False,
                                                transform=ccrs.PlateCarree(),
                                                vmin = 45,
                                                vmax = 105,
                                                )
    
    ax.set_title("")
    
    cbar, ax_cb = add_aligned_colorbar(fig = fig, ax = ax, pc = pc, 
                                        )
    
    cbar.set_label(label = "length of ${V. vulnificus}$ \n "
            "season 2021 [days]",size = 8)
    
    fig, ax = create_map(fig = fig,
                          sp_x = 2, 
                          sp_y = 2, 
                          sp_n = 3,
                          extent = extent,
                          gridticksize = 8)
    add_subplot_char(ax, "c", loc = "lower left")
    ax.ylabel_style = {'size': 8}
    ax.xlabel_style = {'size': 8}
    pc = slope.where(slope != 0).plot(ax = ax,
                                      add_colorbar=False,
                                      transform=ccrs.PlateCarree(),
                                      vmin = -1.5,
                                      vmax = 1.5,
                                      cmap = cm.balance)
    ax.set_title("")
    
    cbar, ax_cb = add_aligned_colorbar(fig = fig, ax = ax, pc = pc)
    
    cbar.set_label(label = "Sen's slope \n[days a$^{-1}$]",size = 8)
    
    fig, ax = create_map(fig = fig,
                          sp_x = 2, 
                          sp_y = 2, 
                          sp_n = 4,
                          extent = extent,
                          left_labels = False,
                          gridticksize = 8)
    
    add_subplot_char(ax, "d", loc = "lower left")
    
    pc = p.where(slope != 0).plot(ax = ax,
                                  add_colorbar=False,
                                  transform=ccrs.PlateCarree(),
                                  cmap = crop_by_percent(cm.oxy, 
                                                          20, 
                                                          which='max', 
                                                          N=None),
                                  vmin = 0,
                                  vmax = 0.2)
    ax.set_title("")
    
    cbar, ax_cb = add_aligned_colorbar(fig = fig, ax = ax, pc = pc, 
                                        label = "p [-]",
                                        ticks = [0, 0.05, 0.1, 0.2])
    cbar.set_label(label = "p [-]",size = 8)
    
    fig.subplots_adjust(wspace = 0.3, hspace = -0.15)
    
    fig.savefig("N:/plots/mann_kendall_sst_17.6_w_lag.jpeg",
                dpi = 300,
                bbox_inches = "tight")
    
    plt.show()
