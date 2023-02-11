# -*- coding: utf-8 -*-
"""
Downloads the "Baltic Sea Biogeochemistry Reanalysis" from CMEMS via OPeNDAP
https://resources.marine.copernicus.eu/product-detail/BALTICSEA_REANALYSIS_BIO_003_012/INFORMATION

Needs PyDAP >= v3.3.0
"""

import xarray as xr
import os
import numpy as np
from tqdm import tqdm
from pydap.client import open_url
from pydap.cas.get_cookies import setup_session

def copernicusmarine_datastore(dataset, username, password):
    """
    Helper function to access the CMEMS store. Modified from CMEMS:
    https://help.marine.copernicus.eu/en/articles/5182598-how-to-consume-the-opendap-api-and-cas-sso-using-python
    """
    cas_url = 'https://cmems-cas.cls.fr/cas/login'
    session = setup_session(cas_url, username, password)
    session.cookies.set("CASTGC", session.cookies.get_dict()['CASTGC'])
    database = ['my', 'nrt']
    url = f'https://{database[0]}.cmems-du.eu/thredds/dodsC/{dataset}'
    try:
        data_store = xr.backends.PydapDataStore(open_url(url, 
                                                         session=session, 
                                                         user_charset='utf-8')) 
    except:
        url = f'https://{database[1]}.cmems-du.eu/thredds/dodsC/{dataset}'
        data_store = xr.backends.PydapDataStore(open_url(url, 
                                                         session=session, 
                                                         user_charset='utf-8'))
    return data_store

os.chdir("N:/data/nutrients/")

username = "" # your CMEMS username
password = "" # your CMEMS password

years = np.arange(2008,2021)

#%% download 4 km Nutrients data since 1993 

dataset_id = 'dataset-reanalysis-scobi-dailymeans'

data_store = copernicusmarine_datastore(dataset_id, 
                                        username = username, 
                                        password = password)

ds = xr.open_dataset(data_store)
ds = ds.sel(latitude = slice(53.6287, 54.95), longitude = slice(9.35, 14.685))

# download files for single years only to avoid runnning into sever errors
for year in tqdm(years):
    ds_year = ds.sel(time = slice("{}-01-01".format(year), 
                                  "{}-12-31".format(year)))
    
    ds_year = ds_year.isel(depth = 0) # select the surface layer
    
    ds_year.to_netcdf("CMEMS_BALTICSEA_REANALYSIS_scobi_{}.nc".format(year))

#%% merge all individual files

files = [x for x in os.listdir() if ".nc" in x]

for i, file in enumerate(files):
    if i == 0:
        ds = xr.open_dataset(file)
        i += 1
    else:
        ds   = xr.concat([ds, xr.open_dataset(file)], dim = "time")
    
ds.to_netcdf("CMEMS_BALTICSEA_REANALYSIS_scobi_all.nc")