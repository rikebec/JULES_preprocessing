# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 17:32:08 2024

@author: Rike Becker

Description: code to create the PTF-input data for JULES from the previously clipped global PFT-ESACCI map 
and reproject the map to the LCC-WRF grid
"""

#%% load libraries
import geopandas as gpd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

#%% 
ptf_path = 'C:\\Users\\uk083025\\Documents\\Imperial_College_London\\Deplete_and_Retreat\\JULES\\JULES_preprocessing\\Output\\netcdf\\pft_clipped_d03_wgs84.nc'
catchments_lcc_path = 'C:\\Users\\uk083025\\Documents\\Imperial_College_London\\Deplete_and_Retreat\\JULES\\JULES_preprocessing\\Output\\shapefiles\\Domain03_catchments_lcc.shp'
catchments_wgs84_path = 'C:\\Users\\uk083025\\Documents\\Imperial_College_London\\Deplete_and_Retreat\\JULES\\JULES_preprocessing\\Input\\shapefiles\\Domain03_small_clip.shp'

#%% set the path to the input data
pft_path = "~/JULES_preprocessing/Output/netcdf/pft_clipped_d03_wgs84.nc"
catchments_lcc_path = "~/JULES_preprocessing/Output/shapefiles/Domain03_catchments_lcc.shp"
catchments_wgs84_path = "~/JULES_preprocessing/Input/shapefiles/Domain03_small_clip.shp"

#%% read the WGS84-catchment data
catchments_wgs84 = gpd.read_file(catchments_wgs84_path, crs="epsg:4326") # read shapefile

#%% read the clipped PFT map
pft = xr.open_dataset(pft_path) 
pft_water = pft.WATER # just for plotting

#%% plot the WGS84 maps to check if catchments are in place
f,ax = plt.subplots()
pft_water.plot(ax=ax)
catchments_wgs84.boundary.plot(ax=ax, color="red")
plt.show()

#%% calculate one PFT variable with fractions of all variables
list(pft.keys())
pft = pft.drop_vars(["spatial_ref", "LAND", "WATER", "WATER_OCEAN"])
water_values = np.unique(np.array(pft_water)) #all values are integers from 1-100
pft_frac = pft/100 #to get fractions

#%% change order of layer (first all veg PFTs, then all non-veg PFTs)
# order of PFTs similar to Harper, A. et al (2016)
# 1. Tropical broadleaf evergreen trees (BET-Tr)
# 2. Temperate broadleaf evergreen trees (BET-Te)
# 3. Broadleaf deciduous trees (BDT)
# 4. Needle-leaf evergreen trees (NET)
# 5. Needle-leaf deciduous trees (NDT)
# 6. C3 Grass
# 7. C4 Grass
# 8. Evergreen Shrubs (ESh)
# 9. Deciduous Shrubs (DSh)
# 10. Urban
# 11. Inland water
# 12. Bare soil
# 13. Land-ice (bearbeitet) 

BARE = np.array(pft.BARE) #13
BUILT = np.array(pft.BUILT) #11
GRASS_MAN = np.array(pft['GRASS-MAN']) #6
GRASS_NAT = np.array(pft['GRASS-NAT']) #5
SHRUBS_BD = np.array(pft['SHRUBS-BD']) #8
SHRUBS_BE = np.array(pft['SHRUBS-BE']) #7
SHRUBS_ND = np.array(pft['SHRUBS-ND']) #10
SHRUBS_NE = np.array(pft['SHRUBS-NE']) #9
WATER_INLAND = np.array(pft.WATER_INLAND) #12
SNOWICE = np.array(pft.SNOWICE) #14
TREES_BD = np.array(pft['TREES-BD']) #2
TREES_BE = np.array(pft['TREES-BE']) #1
TREES_ND = np.array(pft['TREES-ND']) #4
TREES_NE = np.array(pft['TREES-NE']) #3

# join arrays in correct order
layer_order = [TREES_BE, TREES_BD, TREES_NE, TREES_ND, GRASS_NAT, GRASS_MAN, SHRUBS_BE, SHRUBS_BD, SHRUBS_NE, SHRUBS_ND, BUILT, WATER_INLAND, BARE, SNOWICE]

# stack layers
layer_stack = np.stack(layer_order, axis=2)

#%% set SNOWICE == BARE (soil) in pft_param.nml file

#%% create netCDF
# get lat lon values
ds_lat = pft['lat'] 
ds_lon = pft['lon']
# get dimensions
lat_dim = np.array(ds_lat.lat)
lon_dim = np.array(ds_lon.lon)
#get coordinates
lat = np.array(ds_lat)
lon = np.array(ds_lon)

#%%
Nc_img = xr.Dataset(
    coords = {'lon': lon_dim, 'lat': lat_dim, 'type': np.arange(1,len(layer_order)+1)},
    data_vars = {
    'longitude': xr.DataArray(
        data = lon,
        dims = ['lat','lon']
        ),
    'latitude': xr.DataArray(
        data = lat,
        dims = ['lat','lon']
        ),    
    'frac': xr.DataArray(
        data = layer_stack, 
        dims = ['lat', 'lon', 'type']
        ),
    },
   )

#%%
f,ax = plt.subplots()
Nc_img.frac.isel(type=13).plot(ax=ax)
catchments_wgs84.boundary.plot(ax=ax, color="red")
plt.show()

#%% save land cover fraction information as netCDF file to be used in JULES
Nc_img.to_netcdf('C:\\Users\\uk083025\\Documents\\Imperial_College_London\\Deplete_and_Retreat\\JULES\\JULES_preprocessing\\Output\\netcdf\\jules_pft_14_wgs84_dar_d03.nc')
#Nc_img.to_netcdf("~/JULES_preprocesing/Output/netcdf/jules_pft_14_wgs84_dar_d03.nc")


#%% 
#########################################
######## reproject to WRF grid ##########
#########################################

#%% projection:
lcc_crs = {'proj': 'lcc',
           'lat_1': -5,
           'lat_2': -40,
           'lat_0': -22.5000152587891,
           'lon_0': -59.5,
           'x_0': 0,
           'y_0': 0,
           'datum': 'WGS84',
           'units': 'm',
           'no_defs': True}

#%% read the LCC-catchment data
catchments_lcc = gpd.read_file(catchments_lcc_path) 
catchments_lcc.crs = lcc_crs

#%% upscale / resample to WRF grid size


#%% re-project the file to LCC (WRF coordinate system)
pft_calc_proj = pft_calc.to_crs(lcc_crs)
pft_calc.crs

#%% overlay the catchment shapefile (in lcc) to check if they are located where they should be
f,ax = plt.subplots()
pft_water.plot(ax=ax)
catchments_wgs84.boundary.plot(ax=ax, color="k")
plt.show()

#%% save %-map as new netCDF file "jules_pft_2020_d03.nc"


# 1. check if sum of types adds up to 100 or 1
# 2. add coordinates
# 3. reproject
# 4. 
# 5. set up new git repository
# 6. MC fellowship now open - first ideas
# 7. soil maps for jules


# questions for summer school
# 1. use of new PFT map - parameterization information? examples?
# 2. grid - irregular (WRF) grid
# 3. information on data pre-processing
# 4. how to: only simulate single catchments within a larger domain? - set to zero and one?
# 5. how to modify a rose-suite?
# 6. list of files to modify when adopting a rose-suite?

