# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 08:35:41 2024

@author: Rike Becker

Description: script to create latlon-input file for JULES for DAR-domain 3 (Chile + parts of Argentina)

For more information on WRF output coordinate systems see:
    https://fabienmaussion.info/2018/01/06/wrf-projection/
    https://salem.readthedocs.io/en/stable/gis.html
    https://github.com/fmaussion/salem/issues/18

"""
#%% load required libraries
import xarray as xr
import numpy as np
import salem # to be used to read (parse) WRF output files

#%% load WRF output file to get the grid dimensions/number of cells
img_name = "~/JULES_preprocessing/Input/netcdf/wrfout_d03.nc" 

#%% read the image 
ds = salem.open_wrf_dataset(img_name)
#list(ds_salem.keys())
ds.salem.grid

#%%plot landmask
ds.LANDMASK.isel(time=1).salem.quick_map()

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
    'latitude_longitude': xr.DataArray(
        data = -2147483647, # no idea what this number is. 
        attrs = {'grid_mapping_name': 'latitude_longitude', 
                 'longitude_of_prime_meridian': '0.0',
                 'earth_radius': '6371229.0'}
        ),
    },
   )

Nc_img = Nc_img.assign_attrs(pyproj_srs = attribute)
#%% save lat and lon information as netCDF file to be used in JULES
Nc_img.to_netcdf("~/JULES_preprocesing/Output/netcdf/jules_latlong_dar_d03.nc")

#%% plot to check extent
test = Nc_img.longitude
test.plot()
