# -*- coding: utf-8 -*-
"""
Merges the two Vibrio datasets from MV and SH. To do so, the measurements from 
MV are converted into the same units as the measurements from SH, the internet
is scraped for the coordinates of the official bathing sites of MV and finally
the measurements from SH are transformed into the same categorical system as in
MV.
"""

import pandas as pd
import os
import numpy as np
from urllib.request import urlopen
import re
from pyproj import CRS, Transformer
from math import floor, log

os.chdir("N:/data/vibrio/")

#%% import Vibrio data from MV and bring it into the same shape as the SH data

df = pd.read_excel("intermediate/Vibrio_MV_2008-2022.xlsx", 
                   sheet_name = "2008")

for year in np.arange(2009, 2023):

    df_temp = pd.read_excel("intermediate/Vibrio_MV_2008-2022.xlsx", 
                            sheet_name = str(year))
    
    df = pd.concat([df, df_temp])

for char in [".", " ", "l", "L", "/"]:
    df['Konzentration'] = df['Konzentration'].str.replace(char, '', 
                                                          regex = False)

df["quantity"] = np.nan
df.loc[df["Konzentration"] == "negativ", "quantity"] = 0
df.loc[df["Konzentration"] == "positiv", "quantity"] = 1

# Conversion of the MV data from CFU/L to CFU/ml (Faktor 0.001)
df.loc[df["Konzentration"] == "10000", "quantity"] = 10
df.loc[df["Konzentration"] == "100000", "quantity"] = 100
df.loc[df["Konzentration"] == "1000000", "quantity"] = 1000
df.loc[df["Konzentration"] == "10000000", "quantity"] = 10000
df.loc[df["Konzentration"] == "100000000", "quantity"] = 100000

bathing_site_ids = [id for id in pd.unique(df["Badestelle"]) 
                    if type(id) == int]

bathing_sites = []
for i, row in df.iterrows():
    bathing_site = row["Badestelle"]
    if type(bathing_site) == int:
        bathing_site = "DEMV_PR_1_0" + str(bathing_site)
    bathing_sites.append(bathing_site)

df["Badestelle"] = bathing_sites

#%% Scrape the internet for coordinates of official measuring points (as they 
# were not provided by MV, even not upon request) and export the final dataset

def scrape_html(start_str, end_str):
    start_index = html.find(start_str)
    start_index = start_index + len(start_str)
    end_index = html.find(end_str)
    out_str = html[start_index:end_index]
    return out_str

data = []

crs_4326 = CRS("WGS84")
crs_proj = CRS("EPSG:4647")
transformer = Transformer.from_crs(crs_4326, crs_proj)

for bathing_site in bathing_site_ids:
    
    print(bathing_site)

    link = "https://www.regierung-mv.de/Landesregierung/sm/gesundheit/"\
           "Badewasserqualitaet/badewasserkarte/badestelle/"\
           "badegewaesserprofil?gaia.badewasserprofil.id={}"\
           "".format(bathing_site)

    page = urlopen(link)
    html_bytes = page.read()
    html = html_bytes.decode("utf-8")
    
    position = scrape_html(start_str = "<td>Lage (ETRS89)</td><td>", 
                           end_str = "</td></tr></table><!-- setzen wir als h1"
                                     " im Channel")
    
    position = re.findall(r'\d+', position)
    lon = float("{}.{}".format(position[0], position[1]))
    lat = float("{}.{}".format(position[2], position[3]))
    
    name = scrape_html(start_str = "<td>Name des Badegewässers</td><td>", 
                       end_str = "</td></tr><tr><td>Gemeinde</td><td>")
    
    id_number = scrape_html(start_str = "<td>ID Nummer</td><td>", 
                            end_str = "</td></tr><tr><td>Lage (ETRS89)</td>"
                                      "<td>")
    
    county = scrape_html(start_str = "<td>Gemeinde</td><td>", 
                         end_str = "</td></tr><tr>"
                                   "<td>Eigentümer des Gewässers</td><td>")    
    
    east, north = transformer.transform(lat, 
                                        lon)

    data.append({"station_id": id_number,
                 "station_name": name,
                 "typ": np.nan,
                 "kategorie": np.nan,
                 "kreis": county,
                 "east": east,
                 "north": north,
                 "lat": lat,
                 "lon": lon})

# manually add sampling stations which were not at official bathing sites
data.append({"station_id": "HRO17",
             "station_name": "Rostock Hafen LP60S",
             "typ": np.nan,
             "kategorie": np.nan,
             "kreis": "Hansedtadt Rostock",
             "east": 32702252.225625288,
             "north": 6003075.77154479,
             "lat": 54.1359301,
             "lon": 12.0960518})

data.append({"station_id": "HRO18",
             "station_name": "Rostock Stadthafen",
             "typ": np.nan,
             "kategorie": np.nan,
             "kreis": "Hansedtadt Rostock",
             "east": 32704695.143901344,
             "north": 5998424.898305747,
             "lat": 54.0932249,
             "lon": 12.1302371})

data.append({"station_id": "HRO19",
             "station_name": "Rostock-Warnemünde Mittelmole",
             "typ": np.nan,
             "kategorie": np.nan,
             "kreis": "Hansedtadt Rostock",
             "east": 32701646.8120454,
             "north": 6008064.853781942,
             "lat": 54.180943,
             "lon": 12.090136})

data.append({"station_id": "Breitling/Schnater",
             "station_name": "Rostock Breitling am Schnatermann",
             "typ": np.nan,
             "kategorie": np.nan,
             "kreis": "Hansedtadt Rostock",
             "east": 32705034.65590057,
             "north": 6007391.641131102,
             "lat": 54.1735594,
             "lon": 12.1415091})

sites = pd.DataFrame(data)

sites.to_csv("cleaned/Badegewaesser_Ostsee_MV_Koordinaten_cleaned.csv", 
             sep = "\t")

#%% merge coordinates with Vibrio df

vibrio_mv = pd.merge(df, sites, 
                     left_on = "Badestelle", 
                     right_on = "station_id")

vibrio_mv["species"] = "v. vuln."
vibrio_mv["comments"] = np.nan

vibrio_mv = vibrio_mv.rename({#"Badestelle": "station_id",
                              "Nr. ": "old_id",
                              "Datum": "date",
                              "Wassertemperatur": "temperature",
                              "Salinität": "salinity",
                              "Lufttemperatur": "air temp",
                              # "gewässer_name": "station_name",
                              "ostwert": "east"}, axis = 1)

vibrio_mv = vibrio_mv[["date", "station_id", "station_name", "species", 
                       "quantity", "temperature", "air temp", "salinity", 
                       "east", "north", "lat", "lon",
                       "old_id", "comments"
                       ]]

vibrio_mv.to_csv("cleaned/vibrio_data_cleaned_mv.csv", sep = "\t")


#%% merge with data from SH

vibrio_sh = pd.read_csv("cleaned/vibrio_data_cleaned_sh.csv", 
                        sep = "\t",
                        parse_dates = ["date"])

vibrio_df = pd.concat([vibrio_sh, vibrio_mv])
vibrio_df = vibrio_df[vibrio_df["species"] == "v. vuln."]
vibrio_df = vibrio_df.sort_values("date").reset_index(drop = True)

#Transform SH quantities to MV Format
def floor_power_of_10(n):
    if n == 0 or np.isnan(n):
        return n
    else:
        exp = log(n, 10)
        exp = floor(exp)
        return 10**exp

vibrio_df["quantity"] = vibrio_df["quantity"].map(floor_power_of_10)

# drop a duplicate sample from Binz
vibrio_df = vibrio_df.drop(97)

vibrio_df.to_csv("cleaned/vibrio_vuln_cleaned_sh_mv_categorical.csv", 
                 sep = "\t")
