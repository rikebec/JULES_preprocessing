# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:28:06 2024

@author: Rike Becker

Descrpition: script to clip ESA ICC plant functional types from the year 2020 to DAR domains
Requires large memory space to read and process 300m-resolution global map. Hence, run on server or HPC!
source: Harper et al. 2023: https://dx.doi.org/10.5285/26a0f46c95ee4c29b5c650b129aab788
"""
#%% install required libraries
import xarray as xr
import geopandas as gpd
from shapely.geometry import box
import rioxarray as rio #keep even though Spyder might indicate package is not used
#import matplotlib.pyplot as plt
import numpy as np

#%% download the PFT data from https://data.ceda.ac.uk/neodc/esacci/land_cover/data/pft/v2.0.8/ESACCI-LC-L4-PFT-Map-300m-P1Y-2020-v2.0.8.nc

#%% set the paths to the input data
pft_path = "~/JULES_preprocessing/Input/netcdf/ESACCI-LC-L4-PFT-Map-300m-P1Y-2020-v2.0.8.nc"
wrf_path = "~/JULES_preprocessing/Output/netcdf/jules_land_frac_wrf_dar_d04.nc" #to get the extent of the WRF domain
catchments_wgs84_path = "~/JULES_preprocessing/Input/shapefiles/Basins_selected.shp"

#%% read the PFT data
pft = xr.open_dataset(pft_path)
pft = pft.isel(time=0)
pft = pft.drop_vars(["lat_bounds", "lon_bounds", "time_bounds"])
pft.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
pft.rio.write_crs("epsg:4326", inplace=True)

# read the rest of the data
wrf = xr.open_dataset(wrf_path)
catchments_wgs84 = gpd.read_file(catchments_wgs84_path, crs="epsg:4326")

#%% cut region to DAR domains
#get extent of domain
lon_max = np.array(wrf.longitude).max()
lon_min = np.array(wrf.longitude).min()
lat_max = np.array(wrf.latitude).max()
lat_min = np.array(wrf.latitude).min()

# Create a polygon for domain
polygon = box(lon_min, lat_min, lon_max, lat_max)
domain = gpd.GeoDataFrame([1], geometry=[polygon], crs=catchments_wgs84.crs)

#%% clip
pft_clipped = pft.rio.clip(domain.geometry,domain.crs,all_touched=True,drop=True)

#%% save lat and lon information as netCDF file to be used in JULES
pft_clipped.to_netcdf("~/JULES_preprocessing/Output/netcdf/pft_clipped_d04_wgs84.nc")
