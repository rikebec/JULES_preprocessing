# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 09:41:45 2024

@author: Rike Becker

Description: script to create the land fraction input for JULES
"""

#%% load required libraries
import xarray as xr
import numpy as np
import salem

#%% 
# load WRF output file to get the landmask 
img_name = "~/JULES_preprocessing/Input/netcdf/wrfout_d04.nc" 

#%% read the image 
ds = salem.open_wrf_dataset(img_name)
#list(ds.keys())

#%% get landmask data (1= land; 0 = water)
ds_landmask = ds['LANDMASK']
ds_lakemask = ds['LAKEMASK']
land = np.array(ds_landmask.isel(time=1))
lake = np.array(ds_lakemask.isel(time=1))

# make sure that landmask and lakemask add up to 1
land_frac = land+lake
land_frac[land_frac==2]=1

#%% get lat lon values
ds_lat = ds['lat'] 
ds_lon = ds['lon']

# get dimensions
lat_dim = np.array(ds_lat.south_north)
lon_dim = np.array(ds_lon.west_east)

#get coordinates
lat = np.array(ds_lat)
lon = np.array(ds_lon)


#%% get information about the WRF-projection
attribute = ds.pyproj_srs 

#%% create netCDF
Nc_img = xr.Dataset(
    coords = {'lon': lon_dim, 'lat': lat_dim},
    data_vars = {
    'longitude': xr.DataArray(
        data = lon,
        dims = ['lat','lon']
        ),
    'latitude': xr.DataArray(
        data = lat,
        dims = ['lat','lon']
        ),
    'land_frac': xr.DataArray(
        data = land_frac,  
        dims = ['lat','lon']
        )
    },
    )

Nc_img = Nc_img.assign_attrs(pyproj_srs = attribute)
#%% save lat and lon information as netCDF file to be used in JULES
Nc_img.to_netcdf("~/JULES_preprocessing/Output/netcdf/jules_land_frac_lcc_dar_d04.nc")

#%% plot to check extent
test = Nc_img.land_frac
test.plot()

