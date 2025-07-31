#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:56:25 2024

@author: rikebecker

Script to create the slope map (PDM) for JULES .
The data was downloaded and cut to the domains in GEE.
The original DEM used to calculate slope values is the WWF HydroSHEDS void filled DEM (03VFDEM).

"""

#%% load libraries
import xarray as xr
import matplotlib.pyplot as plt
import rioxarray as rio
import rasterio
import numpy as np
import salem
from shapely.geometry import box
import geopandas as gpd
from rasterio.enums import Resampling

#%% Function to upscale the target image to match the reference image resolution
def upscale_to_reference(target_tif, reference_tif, output_tif):
    # Open the reference TIFF file
    with rasterio.open(reference_tif) as ref:
        ref_transform = ref.transform
        ref_crs = ref.crs
        ref_width = ref.width
        ref_height = ref.height
        ref_res = ref.res

    # Open the target TIFF file using rioxarray
    target_ds = rio.open_rasterio(target_tif)
    
    # Reproject and resample the target image to match the reference image resolution
    target_resampled = target_ds.rio.reproject_match(
        rio.open_rasterio(reference_tif),
        resampling=Resampling.bilinear
    )

    # Save the resampled image to a new TIFF file
    target_resampled.rio.to_raster(output_tif)
    print(f"Upscaled image saved to {output_tif}")

#%% 
## ===================================================================== ##
##                         Slope                                         ##
## ===================================================================== ##

#%% set the path to the directory where the input and output data is stored

I_directory = '/Users/rikebecker/Documents/Imperial_College_London/Deplete_and_Retreat/JULES/JULES_preprocessing_IO/Input/earth_engine_exports/d03/'
directory_domain_lcc = '/Users/rikebecker/Documents/Imperial_College_London/Deplete_and_Retreat/JULES/JULES_preprocessing_IO/Output/netcdf/jules_land_frac_wrf_dar_d03.nc'
coords_lcc_path = "/Users/rikebecker/Documents/Imperial_College_London/Deplete_and_Retreat/JULES/JULES_preprocessing_IO/Output/netcdf/jules_land_frac_wrf_dar_d03_clipped.nc"
O_directory = '/Users/rikebecker/Documents/Imperial_College_London/Deplete_and_Retreat/JULES/JULES_preprocessing_IO/Output/netcdf/'

#%% resampling 
file_name = 'slope_d03.tif'
reference_tif = I_directory+'Ksat_d03_0.tif'
target_tif = I_directory+file_name
output_tif = I_directory+'slope_d03_resample.tif'

# perform the upscaling and save resampled image
slope_upscaled = upscale_to_reference(target_tif, reference_tif, output_tif)

# read data
slope_data = rio.open_rasterio(output_tif)
# Remove the extra dimension (squeeze the data)
slope = slope_data.isel(band=0)
slope = slope.drop_vars('band')

slope = slope.rename({'x': 'lon', 'y': 'lat'})
slope.name = 'slope'

# set attributes and missing value (optional)
#slope
slope.attrs['units'] = 'Deg'
slope.attrs['standard_name'] = 'slope'
slope.attrs['long_name'] = 'Slope in Degrees'
#slope.encoding['_FillValue'] = -999.0

#%%
slope.to_netcdf("/Users/rikebecker/Documents/Imperial_College_London/Deplete_and_Retreat/JULES/JULES_preprocessing_IO/Output/netcdf/Slope_test_d03.nc")

#%% clip
slope.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
slope.rio.write_crs("epsg:4326", inplace=True)

# read in the land_frac file to get its coordinate system for later transformation
wrf_land_frac = "/Users/rikebecker/Documents/Imperial_College_London/Deplete_and_Retreat/JULES/JULES_preprocessing_IO/Output/netcdf/jules_land_frac_wrf_dar_d03.nc" 
ds = salem.open_wrf_dataset(wrf_land_frac)
ds.salem.grid # check crs

# get extent of DAR domain 03 (catchments in Chile and Argentina)
lon_max = np.array(ds.longitude).max()
lon_min = np.array(ds.longitude).min()
lat_max = np.array(ds.latitude).max()
lat_min = np.array(ds.latitude).min()

# Create a polygon for domain
polygon = box(lon_min, lat_min, lon_max, lat_max)
domain = gpd.GeoDataFrame([1], geometry=[polygon], crs="epsg:4326")

# plot to check extend of domains
fig, (ax1) = plt.subplots(figsize=(12, 8))
slope.plot(ax=ax1)
domain.boundary.plot(ax=ax1, color="red")
plt.show()

img_clipped = slope.rio.clip(domain.geometry,domain.crs,all_touched=True,drop=True)

img_clipped.to_netcdf("/Users/rikebecker/Documents/Imperial_College_London/Deplete_and_Retreat/JULES/JULES_preprocessing_IO/Output/netcdf/jules_slope_wgs84_dar_d03.nc")

#%% reproject
# read wrf lcc grid with Salem package and reproject and resample the PFT map to the WRF grid
# any file in the wrf grid. Needed to get coordinates.

wrf_grid = salem.open_xr_dataset(coords_lcc_path)
img_reproj = wrf_grid.salem.lookup_transform(img_clipped) #uses "average" method for resampling as default

#%% save land cover fraction information as netCDF file to be used in JULES
output = 'jules_slope_lcc_dar_d03.nc'
O_directory = '/Users/rikebecker/Documents/Imperial_College_London/Deplete_and_Retreat/JULES/JULES_preprocessing_IO/Output/netcdf/'
img_reproj.to_netcdf(O_directory+output)

#%% plot
img_reproj.plot()
