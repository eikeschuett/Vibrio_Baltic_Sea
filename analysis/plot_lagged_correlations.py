# -*- coding: utf-8 -*-
"""
Calculates the Spearman Rank correlation between Vibrio abundances and each
time window of all environmental parameters. Results are plotted and saved in
correlation matrices.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns
import cmocean.cm as cmo
import math
from matplotlib.patches import Rectangle
from tqdm import tqdm

def round_decimals_up(number:float, decimals:int=2):
    """
    Returns a value rounded up to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more")
    elif decimals == 0:
        return math.ceil(number)

    factor = 10 ** decimals
    return math.ceil(number * factor) / factor

from scipy.stats import spearmanr

os.chdir("N:/data/merged_vibrio_env/")

files = [
            'vibrio_mean_sst.csv',
            'vibrio_sst_linreg.csv',
            'vibrio_max_sst.csv',
            'vibrio_mean_sal.csv',
            'vibrio_linreg_sal.csv',
            
            'vibrio_mean_chl.csv',
            'vibrio_linreg_chl.csv',
            'vibrio_mean_o2.csv',
            'vibrio_linreg_o2.csv',
            'vibrio_mean_po4.csv',
            
            'vibrio_linreg_po4.csv',
            'vibrio_mean_no3.csv',
            'vibrio_linreg_no3.csv',
            'vibrio_mean_nh4.csv',
            'vibrio_linreg_nh4.csv',
            
            'vibrio_mean_ASWDIFD_S.csv',
            'vibrio_linreg_ASWDIFD_S.csv',
            'vibrio_mean_T_2M_CL.csv',
            'vibrio_linreg_T_2M_CL.csv',
            'vibrio_mean_tp.csv',
            
            'vibrio_linreg_tp.csv',
            'vibrio_mean_ws_mean.csv',
            'vibrio_linreg_ws_mean.csv',
         ]

#import data 
for file in tqdm(files):
    
    df = pd.read_csv(file, sep = "\t")
    
    # keep only samples from before September 2019
    df["date"] = [pd.to_datetime(x[:10]) for x in df["date"]]
    df = df[df["date"] < "2019-09-01"]

    # remove all stations where the GETM model provides unreasonably small 
    # salinity values (i.e. min of sal < 2. See merge_vibrio_env_data.py for 
    # calculation)
    df = df[~df["station_id"].isin(["DEMV_PR_1_0712",
                                    "DEMV_PR_1_0723",
                                    "DESH_PR_0260",
                                    "DESH_PR_0267",
                                    "DESH_PR_0268",
                                    "HRO18"])]
    
    df = df.drop(['Unnamed: 0.1', 'Unnamed: 0', 'date', 'station_id', 
                  'station_name','species', 'temperature', 'east', 'north', 
                  'lat', 'lon', 'old_id', 'comments', 'air temp', 'salinity'], 
                 axis = 1)

#%%
#calculate lagged spearman correlations and check for significance
    corrs = df.corr(method = "spearman")
    scipy_corrs = spearmanr(df, nan_policy = "omit")
    pval = pd.DataFrame(scipy_corrs[1], 
                        columns = corrs.columns, 
                        index = corrs.index)
    
    corr_mat = np.full([31,31], np.nan)
    pval_mat = np.full([31,31], np.nan)
    
    # corr_mat = np.full([28,28], np.nan)
    # pval_mat = np.full([28,28], np.nan)
    
    
    for column in df.columns[1:]:
        start = int(column.split("_")[-2])
        end = int(column.split("_")[-1])
        
        if file == "vibrio_max_sst.csv":
            start = int((start-30) / 10)
            end = int((end-30) / 10)
        
        corr_mat[start, end] = corrs.loc["quantity", column]
        
        pval_mat[start, end] = pval.loc["quantity", column]
        
    pval_mat = pd.DataFrame(pval_mat)
    pval_str = pval_mat.applymap(lambda x: ''.join(
        ['*' for t in [0.001,0.01,0.05] if x<=t]))
    annot_array = pd.DataFrame(corr_mat).round(2).astype(str) + pval_str
    
    idx = np.argwhere(abs(corr_mat) > np.nanquantile(abs(corr_mat), 0.95))
    
    
    #%%
    #Plot lagged correlations
    
    vmax = round_decimals_up(np.nanmax(abs(corr_mat)), 1)
    vmin = -vmax
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    sns.heatmap(corr_mat, ax = ax, cmap = cmo.balance, vmin = vmin, vmax = vmax,
                cbar_kws={'label': 'Spearman R'},
                annot=annot_array, fmt = "", annot_kws={"size": 3},
                )
    
    if file == "vibrio_max_sst.csv":
        xticklabels = yticklabels = np.arange(30, 310, 10)
        ax.set_xticks(np.arange(0.5, 28.5))
        ax.set_yticks(np.arange(0.5, 28.5))
        ax.set_xticklabels(xticklabels)
        ax.set_yticklabels(yticklabels)
        ax.set_xlim([0, 28])
        ax.set_ylim([28, 0])
    
    for i in idx:
        ax.add_patch(Rectangle((i[1], i[0]),
                               1,1, 
                               fill=False, 
                               edgecolor='w', 
                               lw=1))
        
    the_table = plt.table(cellText=[["*", "p: 0.01-0.05"], 
                                    ["**", "p: 0.001-0.01"],
                                    ["***", "p < 0.001"]],
                      colWidths=[0.05, 0.2],
                      cellLoc = "left",
                      loc='lower left',
                      )
    # remove edge line of the table  
    for key, cell in the_table.get_celld().items():
        cell.set_linewidth(0)
    
    max_idx = np.argwhere(abs(corr_mat) == np.nanmax(abs(corr_mat)))[0]
    
    plt.title("{}; max: start: {} end: {}".format(file, 
                                                  max_idx[0], 
                                                  max_idx[1]))
    
    plt.ylabel("Start date (t - x days)")
    plt.xlabel("End date (t - x days)")
    plt.savefig("N:/plots/corr_matrices/corr_mat_{}.png".format(file[:-4]),
                dpi = 300,
                bbox_inches = "tight")
    plt.show()
    