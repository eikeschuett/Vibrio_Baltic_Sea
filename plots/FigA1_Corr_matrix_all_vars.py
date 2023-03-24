# -*- coding: utf-8 -*-
"""
Plot the correlation matrix between all variables for the appendix (Figure A1).
"""


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr
sns.set_theme(style="white", font_scale = 0.55)

# Import the timelag data to calculate the correlation matrix
fpath = "N:/data/merged_vibrio_env/"
df = pd.read_csv(fpath + "vibrio_data_merged_sh_mv_timelag.csv", 
                 sep = "\t")

# hardcode labels in the correct notation in the order of descending
# absolute r_s with Vibrio qty
labels = ['Vibrio\nqty', 
          'SST\nmean\n0-11', 
          'SAT\nmean\n0-16', 
          'O2\nmean\n0-0',
          'SSS\nmean\n0-6', 
          'SST\n180', 
          'WS\nmean\n2-17',
          'NO3\nmean\n7-25', 
          'Prec\nmean\n6-18', 
          'PO4\nmean\n22-22',
          'SAT\ntrend\n7-29', 
          'Prec\ntrend\n5-8', 
          'Irrad\nmean\n27-28',
          'NO3\ntrend\n0-24', 
          'O2\ntrend\n5-16', 
          'SSS\ntrend\n11-26', 
          'Chl\ntrend\n3-5',
          'WS\ntrend\n4-8', 
          'Irrad\ntrend\n18-28', 
          'PO4\ntrend\n4-9',
          'SST\ntrend\n0-30', 
          'NH4\ntrend\n9-10', 
          'NH4\nmean\n28-28',
          'Chl\nmean\n6-15']
          
order = ['quantity', 'sst_mean_0_11', 'T_2M_CL_mean_0_16', 'o2_mean_0_0',
       'sal_mean_0_6', 'sst_max_180_180', 'ws_mean_mean_2_17',
       'no3_linreg_0_24', 'tp_mean_6_18', 'po4_mean_22_22',
       'T_2M_CL_linreg_7_29', 'tp_linreg_5_8', 'ASWDIFD_S_mean_27_28',
       'no3_mean_7_25', 'o2_linreg_5_16', 'sal_linreg_11_26', 'chl_linreg_3_5',
       'ws_mean_linreg_4_8', 'ASWDIFD_S_linreg_18_28', 'po4_linreg_4_9',
       'sst_linreg_0_30', 'nh4_linreg_9_10', 'nh4_mean_28_28',
       'chl_mean_6_15']

# sort the correlation matrix according the r_s to Vibrio quanity
df = df[order]
rho = df.corr(method = "spearman")
pval = df.corr(method=lambda x, y: spearmanr(x, y)[1]) - np.eye(*rho.shape)
p = pval.applymap(lambda x: ''.join(['*' for t in [.05, .01, .001] if x<=t]))
dat = rho.round(2).astype(str) + p         

# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(rho, dtype=bool))

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(12, 14))

# Generate a custom diverging colormap
cmap = sns.color_palette("vlag", as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
p = sns.heatmap(rho,
                mask=mask,
                cmap=cmap,
                annot=dat,
                # annot_kws={"fontsize":8},
                square = False, 
                linewidths=.5 ,
                vmin = -0.8,
                vmax = 0.8,    
                cbar_kws={"shrink": .5,
                          'label': "$r_s$", 
                          "use_gridspec": False,
                          "location":"top"},
                fmt='')

cbar = ax.collections[0].colorbar
# increase label size
cbar.ax.tick_params(labelsize=12.5)
ax.figure.axes[-1].xaxis.label.set_size(15)

# add a rectangle to highlight the Vibrio column
ax.add_patch(Rectangle((0, 0), 1, 24, fill=False, edgecolor='k', lw=1.5, clip_on=False))

ax.set_xticks(np.arange(0.5, 24.5))
ax.set_yticks(np.arange(0.5, 24.5))
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
ax.set_xlim([0, 24])
ax.set_ylim([24, 0])

plt.xticks(rotation=0)


plt.savefig("N:/plots/corr_matrices/All_vars_Corr.pdf",
            bbox_inches = "tight",
            dpi = 300)