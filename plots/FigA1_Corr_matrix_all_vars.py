# -*- coding: utf-8 -*-
"""
Plot the correlation matrix between all variables for the appendix (Figure A1).
"""


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
sns.set_theme(style="white", font_scale = 0.55)

# Generate a large random dataset
fpath = "N:/data/merged_vibrio_env/"
df = pd.read_csv(fpath + "vibrio_data_merged_sh_mv_timelag.csv", 
                 sep = "\t")

labels = ['Vibrio\nqty',  
           'Chl\nmean\n6-15',
           'Chl\ntrend\n3-5', 
           'O2\nmean\n0-0',
           'O2\ntrend\n5-16', 
           'NO3\nmean\n7-25', 
           'NO3\ntrend\n0-24', 
           'NH4\nmean\n28-28',
           'NH4\ntrend\n9-10', 
           'PO4\nmean\n22-22', # used to be 3\n3
           'PO4\ntrend\n4-9', # instead of 20\n21
           'SST\nmean\n0-11',
           'SST\ntrend\n0-30', # instead of 24\n25
           'SST\n180', #'sst\nmax\n180\n180', 
           'SSS\nmean\n0-6', # sal\nmean\n0\n6
           'SSS\ntrend\n11-26', # sal\ntrend\n11\n26
           'Irrad\nmean\n27-28', #'ASWDIFD\nS\nmean\n27\n28', 
           'Irrad\ntrend\n18-28', #'ASWDIFD\nS\ntrend\n18\n28',
           'SAT\nmean\n0-16', #T\n2M\nCL\nmean\n0\n16', 
           'SAT\ntrend\n7-29', #'T\n2M\nCL\ntrend\n7\n29',
           'Prec\nmean\n6-18', #'tp\nmean\n6\n18',
           'Prec\ntrend\n5-8', #'tp\ntrend\n5\n8', 
           'WS\nmean\n2-17', #'ws\nmean\nmean\n2\n17', 
           'WS\ntrend\n4-8'
          ]

rho = df.corr(method = "spearman")
pval = df.corr(method=lambda x, y: spearmanr(x, y)[1]) - np.eye(*rho.shape)
p = pval.applymap(lambda x: ''.join(['*' for t in [.05, .01, .001] if x<=t]))
dat = rho.round(2).astype(str) + p         

# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(rho, dtype=bool))

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(12, 12))

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
                cbar_kws={"shrink": .5,
                          'label': "$r_s$", 
                          "use_gridspec": False,
                          "location":"top"},
                fmt='')

ax.set_xticks(np.arange(0.5, 24.5))
ax.set_yticks(np.arange(0.5, 24.5))
ax.set_xticklabels(labels)
ax.set_yticklabels(labels)
ax.set_xlim([0, 24])
ax.set_ylim([24, 0])

plt.xticks(rotation=0)


plt.savefig("N:/plots/corr_matrices/All_vars_Corr.jpeg",
            bbox_inches = "tight",
            dpi = 300)