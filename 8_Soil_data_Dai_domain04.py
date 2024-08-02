#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 20:24:23 2024

@author: rikebecker

Description: script to create the soil input netcdf for JULES. The code extracts, stacks, calculates and resamples the required soil information to match the 4x4km WRF grid for Deplete and Retreat
and saves the data as netcdf-file.

Soil Albedo is taken from Houldcroft et al. (2009) https://doi.org/10.1175/2008JHM1021.1

This script uses the dataset of Dai et al. 2019
(https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019MS001784) for soil thermal properties 
(heat capacity 'hcap' and thermal conductivity 'hcon' parameters in JULES)

The soil thermal data was downloaded from: http://globalchange.bnu.edu.cn/research/soil5.jsp
    and pre-processed using the code 'Read_binay_soil_data.py' (see Git repository). 
    This data already contains data for 4 JULES soil layers.

For all other soil parameters the script reads the HiHydroSoil-data and SoilGrids-data, which was clipped to the DAR domain 
and upscaled from 30 m x 30 m to 1 km x 1 km using Google Earth Engine (see Colab notebook in Git repository). 

Source of original data for all soil parameters except for the thermal properties: 
HiHydroSoil v2.0 data on Google Earth Engine: https://gee-community-catalog.org/projects/hihydro_soil/
SoilGrids v2.0 data on Google Earth Engine: https://gee-community-catalog.org/projects/isric/

This data contains information for 6 soil layers according to the following depths:
    (0-0.05 m, 0.05-0.15 m, 0.15-0.3 m, 0.3-0.6 m, 0.6-1.0 m, 1-2 m)
We take the average of layers 1+2 and the average of layers 4+5 to match the 4 JULES layers with the following depths:
JULES layers (0 - 0.1 m, 0.1 - 0.35 m, 0.35 - 1.0 m, and 1.0 - 3.0 m)

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


# Function to safely compute reciprocal of Van Genuchten Alpha value
def safe_reciprocal(array):
    with np.errstate(divide='ignore', invalid='ignore'):
        reciprocal = np.divide(1.0, array)
        reciprocal[~np.isfinite(reciprocal)] = 0  # Replace inf and NaN with zero
    return reciprocal
    
#%% 
## ===================================================================== ##
##                         Soil Albedo                                   ##
## ===================================================================== ##

#%% set the path to the directory where the input and output data is stored
I_directory_1 = '~/JULES_preprocessing/Input/netcdf/'
I_directory_2 = '~/JULES_preprocessing/Input/earth_engine_exports/d04/'
O_directory = '~/JULES_preprocessing/Output/Soil/'

#%% prepare data for clipping to DAR domain

# get extent of DAR domain 03 which was used to clip the google earth engine maps
# this is slightly bigger than the domains

lon_min = -84.5
lat_min = -21.8
lon_max = -62.0
lat_max = -2.8

# Create a polygon for domain
polygon = box(lon_min, lat_min, lon_max, lat_max)
domain = gpd.GeoDataFrame([1], geometry=[polygon], crs="epsg:4326")

#%% read the soil albedo netcdf, clip it to the domain and resample
soil_alb = I_directory_1+'soil_albedo.nc'

soil_albedo = xr.open_dataset(soil_alb)
soil_albedo.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)

# convert coordinates from degrees east (0 to 360) to conventional lat lon values (-180 to 180)
# Extract the longitude values
longitudes = soil_albedo['longitude'].values
# Convert from 0-360 to -180 to 180
longitudes = np.where(longitudes > 180, longitudes - 360, longitudes)
# Update the longitude values in the dataset
soil_albedo = soil_albedo.assign_coords(longitude=longitudes)

soil_albedo.rio.write_crs("epsg:4326", inplace=True)
soil_albedo= soil_albedo.sortby('longitude')
soil_albedo.soil_albedo.plot()

# clip
albedo_clipped = soil_albedo.rio.clip(domain.geometry,domain.crs,all_touched=True,drop=True)
albedo_clipped.soil_albedo.plot()
albedo_clipped.rio.to_raster(O_directory+'soil_albedo_d04.tif')

# downscale to the resolution of 250m to match the SoilGrids and HiHydro data sets
# resampling 
reference_tif = I_directory_2+'Ksat_d04_0.tif'
target_tif = O_directory+'soil_albedo_d04.tif'
output_tif = O_directory+'soil_albedo_d04_resample.tif'

# perform the upscaling and save resampled image
soil_albedo_scaled = upscale_to_reference(target_tif, reference_tif, output_tif)
# open tif for further use
soil_albedo = rio.open_rasterio(output_tif)

soil_albedo = np.array(soil_albedo)
soil_albedo[np.isnan(soil_albedo)] = 0
soil_albedo = np.reshape(soil_albedo, (2160,2505))

#%% 
## ===================================================================== ##
##                         HiHydro                                       ##
## ===================================================================== ##

#%% set the path to the directory where the input and output data is stored
I_directory = '~/JULES_preprocessing/Input/earth_engine_exports/d04/'
directory_domain_lcc = '~/JULES_preprocessing/Output/netcdf/jules_land_frac_lcc_dar_d04.nc'
coords_lcc_path = "~/JULES_preprocessing/Output/netcdf/jules_land_frac_lcc_dar_d04_clipped.nc"
O_directory = '~/JULES_preprocessing/Output/Soil/'

#%% ### sm_crit ###
file_names = ['Sm_crit_d04_0.tif', 'Sm_crit_d04_1.tif', 'Sm_crit_d04_2.tif', 'Sm_crit_d04_3.tif', 'Sm_crit_d04_4.tif', 'Sm_crit_d04_5.tif']

# loop through all files, read and stack
stacked_data = []
for file in file_names:
    with rasterio.open(I_directory+file) as img:
        data = img.read(1)
        result = data*0.0001 # multiply by 0.0001 to get original value (see documentation of data set for more info)
        result = np.maximum(result, 0) # set all negative values to zero
        stacked_data.append(result)
stacked_data = np.stack(stacked_data)

# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_sm_crit = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

#%% ### sm_sat ###
file_names = ['Sm_sat_d04_0.tif', 'Sm_sat_d04_1.tif', 'Sm_sat_d04_2.tif', 'Sm_sat_d04_3.tif', 'Sm_sat_d04_4.tif', 'Sm_sat_d04_5.tif']

# loop through all files, read and stack
stacked_data = []
for file in file_names:
    with rasterio.open(I_directory+file) as img:
        data = img.read(1)
        result = data*0.0001 # multiply by 0.0001 to get original value (see documentation of data set for more info)
        stacked_data.append(result)
stacked_data = np.stack(stacked_data)

# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_sm_sat = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

#%% ### sm_wilt ###
file_names = ['Sm_wilt_d04_0.tif', 'Sm_wilt_d04_1.tif', 'Sm_wilt_d04_2.tif', 'Sm_wilt_d04_3.tif', 'Sm_wilt_d04_4.tif', 'Sm_wilt_d04_5.tif']

# loop through all files, read and stack
stacked_data = []
for file in file_names:
    with rio.open_rasterio(I_directory+file) as img:
        data = img[0]
        result = data*0.0001 # multiply by 0.0001 to get original value (see documentation of data set for more info)
        stacked_data.append(result)
stacked_data = np.stack(stacked_data)

# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_sm_wilt = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

#%% ### satcon ###
file_names = ['Ksat_d04_0.tif', 'Ksat_d04_1.tif', 'Ksat_d04_2.tif', 'Ksat_d04_3.tif', 'Ksat_d04_4.tif', 'Ksat_d04_5.tif']

# loop through all files, read and stack
stacked_data = []
for file in file_names:
    with rio.open_rasterio(I_directory+file) as img:
        data = img[0]
        result = data*0.0001 # multiply by 0.0001 to get original value (see documentation of data set for more info)
        result =  (result * 10) / (24 * 60 * 60)   # convert from cm/d to kg/m2/s
        stacked_data.append(result)
stacked_data = np.stack(stacked_data)

# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_ksat = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

#%% ### alpha ### in 1/cm
file_names = ['Alpha_d04_0.tif', 'Alpha_d04_1.tif', 'Alpha_d04_2.tif', 'Alpha_d04_3.tif', 'Alpha_d04_4.tif', 'Alpha_d04_5.tif']

# loop through all files, read and stack
stacked_data = []
for file in file_names:
    with rio.open_rasterio(I_directory+file) as img:
        data = img[0]
        result = data*0.0001 # multiply by 0.0001 to get original value (see documentation of data set for more info)
        stacked_data.append(result)
stacked_data = np.stack(stacked_data)

# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_alpha = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

#%% 
## ===================================================================== ##
##                         SoilGrids                                     ##
## ===================================================================== ##

#%% ### SAND #### in g/kg (divide by 10 to get g/100g) - conversion is done at later step
# resampling 
file_name = 'Sand_d04.tif'
reference_tif = I_directory+'Ksat_d04_0.tif'
target_tif = I_directory+file_name
output_tif = I_directory+'sand_d04_resample.tif'

sand = rio.open_rasterio(I_directory+file_name)
sand_0 = np.array(sand[0])

# perform the upscaling and save resampled image
sand_upscaled = upscale_to_reference(target_tif, reference_tif, output_tif)

# Perform the averaging to get 4 soil layers
stacked_data = rio.open_rasterio(output_tif)
# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_sand = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

# set all fill values of -32768 to Nan
stacked_data_sand[stacked_data_sand == -32768] = np.nan

#%% ### SILT ### in g/kg (divide by 10 to get g/100g) - conversion is done at later step
file_name = 'Silt_d04.tif'
reference_tif = I_directory+'Ksat_d04_0.tif'
target_tif = I_directory+file_name
output_tif = I_directory+'silt_d04_resample.tif'

silt = rio.open_rasterio(I_directory+file_name)
silt_0 = np.array(silt[0])

# perform the upscaling and save resampled image
silt_upscaled = upscale_to_reference(target_tif, reference_tif, output_tif)

# Perform the averaging to get 4 soil layers
stacked_data = rio.open_rasterio(output_tif)
# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_silt = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

# set all fill values of -32768 to Nan
stacked_data_silt[stacked_data_silt == -32768] = np.nan

#%% ### CLAY ### in g/kg (divide by 10 to get g/100g) - conversion is done at later step
file_name = 'Clay_d04.tif'
reference_tif = I_directory+'Ksat_d04_0.tif'
target_tif = I_directory+file_name
output_tif = I_directory+'clay_d04_resample.tif'

clay = rio.open_rasterio(I_directory+file_name)
clay_0 = np.array(clay[0])

# perform the upscaling and save resampled image
clay_upscaled = upscale_to_reference(target_tif, reference_tif, output_tif)

# Perform the averaging to get 4 soil layers
stacked_data = rio.open_rasterio(output_tif)
# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_clay = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

# set all fill values of -32768 to Nan
stacked_data_clay[stacked_data_clay == -32768] = np.nan

#%% ### Soil Organic Carbon ### in dg/kg (multiply by 10 to get conventional units in g/kg) - done in PTF
file_name = 'soc_d04.tif'
reference_tif = I_directory+'Ksat_d04_0.tif'
target_tif = I_directory+file_name
output_tif = I_directory+'soc_d04_resample.tif'

soc = rio.open_rasterio(I_directory+file_name)

# perform the upscaling and save resampled image
soc_upscaled = upscale_to_reference(target_tif, reference_tif, output_tif)

# Perform the averaging to get 4 soil layers
stacked_data = rio.open_rasterio(output_tif)
# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_soc = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

# set all fill values of -32768 to Nan
stacked_data_soc[stacked_data_soc == -32768] = np.nan


#%% ### Cation Exchange Capacity ### in mmol/kg (devide by 10 to get conventional units in cmol/kg) - not done in PTF ???
file_name = 'cec_d04.tif'
reference_tif = I_directory+'Ksat_d04_0.tif'
target_tif = I_directory+file_name
output_tif = I_directory+'cec_d04_resample.tif'

cec = rio.open_rasterio(I_directory+file_name)

# perform the upscaling and save resampled image
cec_upscaled = upscale_to_reference(target_tif, reference_tif, output_tif)

# Perform the averaging to get 4 soil layers
stacked_data = rio.open_rasterio(output_tif)
# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_cec = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

# set all fill values of -32768 to Nan
stacked_data_cec[stacked_data_cec == -32768] = np.nan

#%% ### pH ### in pH x 10 (devide by 10 to get conventional units in pH) 
file_name = 'pH_d04.tif'
reference_tif = I_directory+'Ksat_d04_0.tif'
target_tif = I_directory+file_name
output_tif = I_directory+'pH_d04_resample.tif'

pH = rio.open_rasterio(I_directory+file_name)

# perform the upscaling and save resampled image
pH_upscaled = upscale_to_reference(target_tif, reference_tif, output_tif)

# Perform the averaging to get 4 soil layers
stacked_data = rio.open_rasterio(output_tif)
# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_pH = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

# set all fill values of -32768 to Nan
stacked_data_pH[stacked_data_pH == -32768] = np.nan

#%% ### Bulk Density ### in cg/cm³ (devide by 100 to get conventional units in kg/dm³) - done in PTF
file_name = 'bd_d04.tif'
reference_tif = I_directory+'Ksat_d04_0.tif'
target_tif = I_directory+file_name
output_tif = I_directory+'bd_d04_resample.tif'

bd = rio.open_rasterio(I_directory+file_name)

# perform the upscaling and save resampled image
bd_upscaled = upscale_to_reference(target_tif, reference_tif, output_tif)

# Perform the averaging to get 4 soil layers
stacked_data = rio.open_rasterio(output_tif)
# Compute the mean for layers 1-2 and 4-5
average_subset1 = np.mean([stacked_data[0],stacked_data[1]], axis=0)
average_subset2 = np.mean([stacked_data[3],stacked_data[4]], axis=0)
stacked_data_bd = np.stack([average_subset1, stacked_data[2], average_subset2, stacked_data[5]])

# set all fill values of -32768 to Nan
stacked_data_bd[stacked_data_bd == -32768] = np.nan

#%% ### hcap, hcon and brooks-corey "b" coefficient ###
# equations taken from Zed Zulkafli's PhD Thesis

# Units: J m-3 K-1
cc=2373000.00
cs=2133000.00
csi=2133000.00
lambda_air=0.025
lambda_clay=1.16025
lambda_sand=1.57025
lambda_silt=1.57025

clay = stacked_data_clay*0.001 # unit in percentage
sand = stacked_data_sand*0.001
silt = stacked_data_silt*0.001

soil_organic_carbon = stacked_data_soc*0.1
cation_exchange_capacity = stacked_data_cec*0.1
ph_index = stacked_data_pH*0.1
bulk_density = stacked_data_bd*0.01
Ksat = stacked_data_ksat #conversion already done in code above
alpha = stacked_data_alpha 

#theta_sat = 0.01 * (81.799 + (0.099 * clay)
#        - (31.42 * bulk_density) + (0.018 * cation_exchange_capacity)
#        + (0.451 * ph_index / 10) - (0.0005 * sand * clay))

theta_sat = 0.505 - (0.037 * clay) - (0.142 * sand)

# heat capacity (in J m-3 K-1)
hcap = (1 - theta_sat) * (clay * cc + sand * cs + silt * csi)

# dry thermal conductivity (in W m-1 K-1 )
hcon = (lambda_air ** theta_sat) * (lambda_clay ** ((1 - theta_sat) * clay)) * (lambda_sand ** ((1 - theta_sat) * sand)) * (lambda_silt ** ((1 - theta_sat) * silt))   

b = (3.10 + 15.7 * clay - 0.3 * sand) #brooks-corey

# sathh - soil matric suction at saturation (Van Genuchten) in m
#alpha = 9.80665 * np.exp((-2.294 - (3.526 * silt)
#      + (2.440 * (soil_organic_carbon)) - (0.076 * cation_exchange_capacity)
#      - (11.331 * ph_index) + (0.019 * silt * silt)))

sathh = safe_reciprocal(alpha)*0.01 #1/alpha and conversion from cm to m

#%% read one image to get spatial ref and dimensions
dim = rio.open_rasterio(I_directory+file_names[0])
# get coordinates
lon_coord = np.array(dim.x)
lat_coord = np.array(dim.y)

#%% create netCDF 
Nc_img = xr.Dataset(
    coords = {'lon': lon_coord, 'lat': lat_coord, 'soil': np.arange(1,len(stacked_data_sm_crit)+1)},
    data_vars = {   
    'albsoil': xr.DataArray(
        data = soil_albedo, 
        dims = ['lat', 'lon']
        ),      
    'b': xr.DataArray(
        data = b, 
        dims = ['soil', 'lat', 'lon']
        ),    
    'hcap': xr.DataArray(
        data = hcap, 
        dims = ['soil', 'lat', 'lon']
        ),
    'hcon': xr.DataArray(
        data = hcon, 
        dims = ['soil', 'lat', 'lon']
        ),
    'satcon': xr.DataArray(
        data = Ksat, 
        dims = ['soil', 'lat', 'lon']
        ),
    'sathh': xr.DataArray(
        data = sathh, 
        dims = ['soil', 'lat', 'lon']
        ),    
    'sm_crit': xr.DataArray(
        data = stacked_data_sm_crit, 
        dims = ['soil', 'lat', 'lon']
        ),
    'sm_sat': xr.DataArray(
        data = stacked_data_sm_sat, 
        dims = ['soil', 'lat', 'lon']
        ),
    'sm_wilt': xr.DataArray(
        data = stacked_data_sm_wilt, 
        dims = ['soil', 'lat', 'lon']
        ),
    },
   )

# set attributes and missing value (optional)
#Soil Albedo
Nc_img.albsoil.attrs['units'] = '1'
Nc_img.albsoil.attrs['standard_name'] = 'albsoil'
Nc_img.albsoil.attrs['long_name'] = 'Soil_Albedo'
Nc_img.albsoil.encoding['_FillValue'] = -999.0
#b
Nc_img.b.attrs['units'] = '1'
Nc_img.b.attrs['standard_name'] = 'b'
Nc_img.b.attrs['long_name'] = 'exponent in soil hydraulic characteristics'
Nc_img.b.encoding['_FillValue'] = -999.0
#hcap
Nc_img.hcap.attrs['units'] = 'J m-3 K-1'
Nc_img.hcap.attrs['standard_name'] = 'hcap'
Nc_img.hcap.attrs['long_name'] = 'dry heat capacity'
Nc_img.hcap.encoding['_FillValue'] = -999.0
#hcon
Nc_img.hcon.attrs['units'] = 'W m-1 K-1'
Nc_img.hcon.attrs['standard_name'] = 'hcon'
Nc_img.hcon.attrs['long_name'] = 'dry thermal conductivity'
Nc_img.hcon.encoding['_FillValue'] = -999.0
#satcon
Nc_img.satcon.attrs['units'] = 'kg m-2 s-1'
Nc_img.satcon.attrs['standard_name'] = 'satcon'
Nc_img.satcon.attrs['long_name'] = 'hydraulic conductivity at saturation'
Nc_img.satcon.encoding['_FillValue'] = -999.0
#sathh
Nc_img.sathh.attrs['units'] = 'm'
Nc_img.sathh.attrs['standard_name'] = 'sathh'
Nc_img.sathh.attrs['long_name'] = 'absolute value of the soil matric suction at saturation'
Nc_img.sathh.encoding['_FillValue'] = -999.0
#sm_crit
Nc_img.sm_crit.attrs['units'] = 'm3 m-3'
Nc_img.sm_crit.attrs['standard_name'] = 'sm_crit'
Nc_img.sm_crit.attrs['long_name'] = 'volumetric soil moisture content at critical point'
Nc_img.sm_crit.encoding['_FillValue'] = -999.0
#sm_sat
Nc_img.sm_sat.attrs['units'] = 'm3 m-3'
Nc_img.sm_sat.attrs['standard_name'] = 'sm_sat'
Nc_img.sm_sat.attrs['long_name'] = 'volumetric soil moisture content at saturation'
Nc_img.sm_sat.encoding['_FillValue'] = -999.0
#sm_wilt
Nc_img.sm_wilt.attrs['units'] = 'm3 m-3'
Nc_img.sm_wilt.attrs['standard_name'] = 'sm_wilt'
Nc_img.sm_wilt.attrs['long_name'] = 'volumetric soil moisture content at wilting point'
Nc_img.sm_wilt.encoding['_FillValue'] = -999.0

Nc_img.to_netcdf("~/JULES_preprocessing/Output/netcdf/Soil_test_d04.nc")

#%% clip
Nc_img.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
Nc_img.rio.write_crs("epsg:4326", inplace=True)

# read in the land_frac file to get its coordinate system for later transformation
wrf_land_frac = "~/JULES_preprocessing/Output/netcdf/jules_land_frac_lcc_dar_d04.nc" 
ds = salem.open_wrf_dataset(wrf_land_frac)
ds.salem.grid # check crs

# get extent of DAR domain 04 (catchments in Peru and Bolivia)
lon_max = np.array(ds.longitude).max()
lon_min = np.array(ds.longitude).min()
lat_max = np.array(ds.latitude).max()
lat_min = np.array(ds.latitude).min()

# Create a polygon for domain
polygon = box(lon_min, lat_min, lon_max, lat_max)
domain = gpd.GeoDataFrame([1], geometry=[polygon], crs="epsg:4326")

# plot to check extend of domains
fig, (ax1) = plt.subplots(figsize=(12, 8))
Nc_img.sm_wilt[0].plot(ax=ax1)
domain.boundary.plot(ax=ax1, color="red")
plt.show()

img_clipped = Nc_img.rio.clip(domain.geometry,domain.crs,all_touched=True,drop=True)

img_clipped.to_netcdf("~/JULES_preprocessing/Output/netcdf/jules_soil_wgs84_dar_d04.nc")

#%% reproject
# read wrf lcc grid with Salem package and reproject and resample the PFT map to the WRF grid
# any file in the wrf grid. Needed to get coordinates.

wrf_grid = salem.open_xr_dataset(coords_lcc_path)
img_reproj = wrf_grid.salem.lookup_transform(img_clipped) #uses "average" method for resampling as default

#%% save as netCDF file to be used in JULES
output = 'jules_soil_lcc_dar_d04.nc'
O_directory = '~/JULES_preprocessing/Output/netcdf/'
img_reproj.to_netcdf(O_directory+output)

#%% plot
test = img_reproj.sm_crit[1]
test.plot()

#%% 
## =================================================================================== ##
##        change soil thermal parameters to the ones from Dai et al. (2019)            ##
## =================================================================================== ##

#%% set the path to the directory where the input and output data is stored
I_directory = '~/JULES_preprocessing/Output/Soil/'
O_directory = '~/JULES_preprocessing/Output/netcdf/'

#%% hcap - dry heat capacity
hcap_dai = I_directory+'csol_d04/hcap_d04_lcc.nc'
hcap_dai = xr.open_dataset(hcap_dai)
hcap_dai = np.array(hcap_dai.hcap)

# plot difference of two data sets. Note!: Difference is large over glacierized areas!!!
soil = xr.open_dataset('~/JULES_preprocessing/Output/netcdf/jules_soil_lcc_dar_d04.nc')
hcap = np.array(soil.hcap)

test = hcap[1] - hcap_dai[1]
plt.imshow(test, interpolation="nearest", origin="upper")
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()

#%% hcon - thermal conductivity
hcon_dai = I_directory+'tkdry_d04/hcon_d04_lcc.nc'
hcon_dai = xr.open_dataset(hcon_dai)
hcon_dai = np.array(hcon_dai.hcon)

# plot difference of two data sets. Note!: Difference is large over glacierized areas!!!
hcon = np.array(soil.hcon)

test = hcon[1] - hcon_dai[1]
plt.imshow(test, interpolation="nearest", origin="upper")
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()

# apply filtering to get rid of irregularities in the image
# select by difference
# set to NAN
# interpolate NAN values

#%% change hcon and hcap to Dai et al. instead of PTF derived params
soil['hcap'].values = hcap_dai
soil['hcon'].values = hcon_dai

#%% save as netCDF file to be used in JULES
output = 'jules_soil_lcc_dar_d04_Dai.nc'
soil.to_netcdf(O_directory+output)
