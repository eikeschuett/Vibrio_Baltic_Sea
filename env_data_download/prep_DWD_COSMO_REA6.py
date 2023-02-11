# -*- coding: utf-8 -*-
"""
Downloads COSMO REA6 data from the German Weather Service (DWD).
Refer to https://opendata.dwd.de/climate_environment/REA/ for more 
information on the reanalysis.
"""

import numpy as np
import os
import urllib.request
from pathlib import Path

os.chdir("N:/data/meteo_data/COSMO_REA6/raw")

months = np.arange(1,13)
years = np.arange(2008, 2020)

# define all the variables for which data shall be downloaded
varnames = ["T_2M", "ASWDIFD_S","U_10M","V_10M", "TOT_PRECIP"] 

url_base = "https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/daily/2D/"

for var in varnames[4:]:
    url = "{}/{}/".format(url_base, var)
    for year in years:
        for month in months:
            if year == 2019 and month > 8:
                continue
            if var == "TOT_PRECIP":
                fname = "{var}.2D.{year}{month:02d}.DaySum.grb".format(
                                                                var = var,
                                                                year = year,
                                                                month = month)
            else:
                fname = "{var}.2D.{year}{month:02d}.DayMean.grb".format(
                                                                var = var,
                                                                year = year,
                                                                month = month)
            
            Path("").mkdir(parents=True, exist_ok=True)
            if not os.path.exists(var):
                os.makedirs(var)
                
            urllib.request.urlretrieve(url+fname, var + "/" + fname)
