# -*- coding: utf-8 -*-
"""
Creates the plot of the moving window system (Figure 2 at the moment).
"""

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as patches
import numpy as np

bbox = dict(facecolor='w', alpha=0.8, edgecolor = "None")
avg_offset = 0.075

with xr.open_dataset("N:/data/GETM/WB200m.suf.2008-2021.nc") as ds:
    data = ds.sel(lat = 54.348002,
                  lon = 10.165014,
                  time = pd.date_range("2018-07-01", periods=35),
                  method = "nearest")["temp"]
    
    # replace certain very high values
    data[23:26] = [18.4, 19.8, 19.9]
    
    #%%
    xticks = ["2018-07-15", 
              "2018-07-20", "2018-07-25", "2018-07-30"]
    
    fig = plt.figure(figsize = (6.5,3))
    
    ax1 = fig.add_subplot(1,1,1)
    data.plot(ax = ax1)
    date = pd.to_datetime("2018-07-30")
    date2 = pd.to_datetime("2018-07-30T10:00:00")
        
    # ax1.text(date2, 18.2, 'Vibrio measurement', 
    #         {'ha': 'center', 'va': 'center', 'bbox': None},
    #         rotation=90,
    #         fontsize = 6)
    
    
    td_off_str = pd.Timedelta(hours = 24)
    date_end = pd.to_datetime("2018-07-30")
    td_offset = pd.Timedelta(hours = 0)
    for i, td_int in enumerate(np.arange(-5,-20, -5)):
        td = pd.Timedelta(days = td_int)
        date_start = date_end + td
        avg = float(data.sel(time = pd.date_range(date_start, 
                                                  periods = abs(td_int))).mean())
        
        
        rect = patches.Rectangle((date_start-td_offset, avg-0.05), 
                                  abs(td), 0.1,
                                  linewidth=1, 
                                  edgecolor='k', 
                                  facecolor='lightgrey')
        # Add the patch to the Axes
        ax1.add_patch(rect)
        ax1.text(date_start-pd.Timedelta(days = td_int-0.5), 
                 avg-avg_offset, 
                 r"$\overline{SST}_{t" + u"\mathrm{\u2212}" + str(i+1) + ",t" + 
                 u"\mathrm{\u2212}" + "0}$",
                 bbox = bbox)
        
    date_end = pd.to_datetime("2018-07-25")
    
    for i, td_int in enumerate(np.arange(-5,-15, -5)):
        td = pd.Timedelta(days = td_int)
        date_start = date_end + td
        avg = float(data.sel(time = pd.date_range(date_start, 
                                                  periods = abs(td_int))).mean())
        
        
        rect = patches.Rectangle((date_start-td_offset, avg-0.05), 
                                  abs(td), 0.1, 
                                  linewidth=1, 
                                  edgecolor='k', 
                                  facecolor='lightgrey')
        # Add the patch to the Axes
        ax1.add_patch(rect)
        ax1.text(date_start-pd.Timedelta(days = td_int-0.5), 
                 avg-avg_offset, 
                 r"$\overline{SST}_{t" + u"\mathrm{\u2212}" + str(i+2) + ",t" + 
                 u"\mathrm{\u2212}" +"1}$",
                 # bbox = bbox,
                 )
        
    date_end = pd.to_datetime("2018-07-20")
        
    for i, td_int in enumerate(np.arange(-5,-10, -5)):
        td = pd.Timedelta(days = td_int)
        date_start = date_end + td
        avg = float(data.sel(time = pd.date_range(date_start, 
                                                  periods = abs(td_int))).mean())
        
        rect = patches.Rectangle((date_start-td_offset, avg-0.05), 
                                  abs(td), 0.1, 
                                  linewidth=1, 
                                  edgecolor='k', 
                                  facecolor='lightgrey')
        # Add the patch to the Axes
        ax1.add_patch(rect)  
        ax1.text(date_start-pd.Timedelta(days = td_int-0.5), 
                 avg-avg_offset, 
                 r"$\overline{SST}_{t" + u"\mathrm{\u2212}" + str(i+3) + ",t" + 
                 u"\mathrm{\u2212}" + "2}$",
                 bbox = bbox)
    
    ax1.axvline(x=date, c = "Gray")
        
    ax1.set_xticks(xticks)
    ax1.set_xticklabels([])
    ax1.set_yticks([17, 18, 19, 20])
    ax1.set_xlabel("time")
    ax1.set_ylabel("SST [°C]")
    ax1.set_xlim(pd.to_datetime("2018-07-14"), pd.to_datetime("2018-08-03"))
    ax1.set_ylim([17.7, 20.])
    ax1.set_title("")
    
    ax1.set_xticks(xticks)
    ax1.set_xticklabels([u"t\u22123", u"t\u22122", u"t\u22121", u"t\u22120"], 
                        rotation = "horizontal",
                        ha = "center")
    
    plt.savefig("N:/plots/variable_moving_window.pdf",
                dpi = 300,
                bbox_inches = "tight")
    
    plt.show()
    
    