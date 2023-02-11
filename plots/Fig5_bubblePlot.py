# -*- coding: utf-8 -*-8
"""
Creates the bubble plot of Vibrio quanities over SST and SSS (Figure 5).
"""

import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt

# Plot Bubble plot

plt.figure(figsize = (7,4))
data = pd.read_csv('N:/data/merged_vibrio_env/vibrio_data_merged_sh_mv.csv', 
                   parse_dates = ["date"], 
                   sep = "\t")

data["Month"] = data["date"].dt.month
data1 = data.loc[data["quantity"] != 0]
data2 = data.loc[data["quantity"] == 0]
##Apply logarithmic function to the vibrio quantity 
data1['Quantity log10'] =  np.log10(data1['quantity']) 

sns.set(rc={"figure.dpi":600, 'savefig.dpi':400}) 
sns.set_theme(style="ticks")
color1 = ['#66c2a5','#fc8d62','#8da0cb',
          '#e78ac3','#a6d854','#ffd92f','#e5c494']
color = color1[2:-1]

g = sns.scatterplot(data = data1, 
                    x = "sal_mean_0_6", 
                    y = "sst_mean_0_11", 
                    size = "Quantity log10",
                    sizes = (10,250), 
                    hue = "Month", 
                    palette = sns.color_palette(color),
                    alpha=0.75)

x = sns.scatterplot(data = data2, 
                    x = "sal_mean_0_6", 
                    y = "sst_mean_0_11",
                    hue ="Month",
                    marker = "x",
                    palette = sns.color_palette(color1) )

g.set(ylabel = "$\overline{SST}_{0-11} \ $[°C]", 
      xlabel = "$\overline{SSS}_{0-6} \ $ [PSU]")

h,l = g.get_legend_handles_labels()
lalt = ['Quantity', '0.1','1','10','100','10000','100000']
h1,l1 = g.get_legend_handles_labels()


month_l = ["Apr","Mar","Jun","Jul","Aug","Sep","Oct"]

h2 = np.append(np.append(np.append(h1[0], np.append(h1[0],h[6:12])),
                         h1[0]), h[12:19])

# for months as short names
l2 = np.append(np.append( np.append( np.append(lalt[0],
                                               "[CFU $ml^{-1}$]") ,
                                    lalt[1:]),l1[0]), month_l) 

legend = g.legend(h2,l2, 
                  loc = "lower right",
                  frameon=False,
                  ncol = 2, 
                  prop={'size': 8})

plt.savefig("N:/plots/bubble_SST_Sal.jpeg", 
            bbox_inches = "tight",
            dpi = 300)





















