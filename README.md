[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]
# First steps towards a near real-time modelling system of Vibrio vulnificus in the Baltic Sea
This repository contains the code for the manuscript "First steps towards a near real-time modelling system of Vibrio vulnificus in the Baltic Sea", submitted to the International Journal of Environmental Research and Public Health. In this research, we compiled a dataset with more than 600 samples of Vibrio vulnificus and then use data or models and reanalysis to derive the environmental network between several hydrodynamic, meteorological and biogeochemical parameters and V. vulnificus in the south-western Baltic Sea. We also analyse changes of the season length between 1995 and 2021.

Except a manual cleaning of the original Vibrio datasets provided by the health agencies, the data analysis was conducted in Python. The code in this repo can be used to reproduce the results and all figures.

# Workflow:
1. Download and prepare all environmental data
    1. To download the CMEMS data through openDAP run ``env_data_download/prep_CMEMS_nutrients.py``
    1. To download COSMO-REA6 data from the FTP server of the DWD run ``env_data_download/prep_DWD_COSMO_REA6.py``
    1. GETM Data (SST and SSS) is available [on the thredds server of the IOW](https://thredds-iow.io-warnemuende.de/thredds/catalogs/regions/baltic/regions/catalog_WB200m.SST.html)
        1. To validate its SST and SSS products product, run ``validation_of_inputs/validation_GETM_SST_SSS.py``
1. Prepare data of the Vibrio measurements
    1. to merge intermediate data from MV with data from SH run ``data_preparation/prep_mv_data_and_merge_w_sh.py``
    1. 
1. Merge Vibrio data and environmental data and derive the optimal time lag
    1. Apply the flixble moving window of flexible size to calculate means and trends within each window on the GETM data (SST and SSS) and CMEMS nutrients data: run ``data_preparation/merge_vibrio_env_data.py``
    1. To do so for the meteorological data (DWD COSMO REA 6), run ``data_preparation/merge_vibrio_meteo_data.py``
    1. To plot the lagged correlations in correlation matrices, run ``analysis/plot_lagged_correlations.py``. This creates correlation matrices by default in ``N:/plots/corr_matrices/``
    1. Drop stations with unreliable SSS, select optimal time lags based on the correlation matrices, add them into/change them in ``data_preparation/compile_vibrio_and_env_timelag.py`` and run this script. Output CSV files will by default be written to ``N:/data/merged_vibrio_env/``
1. Run the [St. Nicolas House Analysis](https://github.com/thake93/snha4py) with ``analysis/run_snha_for_Vibrios.py``
1. Analyse the SST time series to calculate changes of the Vibrio season length by running ``analysis/Vibrio_season.py``
1. Figures 3, 6, 7 and A2 are plotted from within the functions above. Additional plots were created with code in ``plots/``.

# License


This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
