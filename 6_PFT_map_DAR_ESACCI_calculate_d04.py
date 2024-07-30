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
import salem

#%% set the path to the input data
# previously clipped global PFT map (300m resolution in wgs84 projection)
pft_path = "~/JULES_preprocessing/Output/netcdf/pft_clipped_d04_wgs84.nc"
# shapefile of DAR catchments in the WRF-LCC grid (optional - only for plotting and testing)
catchments_lcc_path = "~/JULES_preprocessing/Output/shapefiles/Domain04_catchments_lcc.shp"
# shapefiles of DAR catchments in WGS84 grid (optional - only for plotting and testing)
catchments_wgs84_path = "~/JULES_preprocessing/Input/shapefiles/Domain04_small_clip.shp"
# any file in the wrf grid. Needed to get coordinates.
coords_lcc_path = "~/JULES_preprocessing/Output/netcdf/jules_land_frac_wrf_dar_d04_clipped.nc"

#%% read the clipped PFT map in WGS84 coordinate system
pft = salem.open_xr_dataset(pft_path) 
pft_water = pft.WATER # optional - just for plotting

#%% (optional) plot the WGS84 maps to check if catchments are in place
# read the WGS84-catchment data
catchments_wgs84 = gpd.read_file(catchments_wgs84_path, crs="epsg:4326") # read shapefile
# plot
f,ax = plt.subplots()
pft_water.plot(ax=ax)
catchments_wgs84.boundary.plot(ax=ax, color="red")
plt.show()

#%% check values and delete the ones which are not needed
list(pft.keys())
pft = pft.drop_vars(["spatial_ref", "LAND", "WATER", "WATER_OCEAN"])
water_values = np.unique(np.array(pft_water)) #all values are integers from 1-100
water_values
pft = pft.reset_coords(['time'], drop=True)

#%% read wrf lcc grid with Salem package and reproject and resample the PFT map to the WRF grid
wrf_grid = salem.open_xr_dataset(coords_lcc_path)
pft_wrf_reproj = wrf_grid.salem.lookup_transform(pft) #uses "average" method for resampling as default

#%% get fractions
pft_frac = pft_wrf_reproj/100 #to get fractions. This increases the data size significantly. Compress when saving to netcdf!
water_values = np.unique(np.array(pft_frac.WATER_INLAND)) # check fractions
water_values
# round to one decimal
pft_frac = pft_frac.round(decimals=1)

#%% (optional) assign crs (lcc) to catchment shapefile and plot the LCC maps to check if catchments are in place
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
catchments_lcc = gpd.read_file(catchments_lcc_path) 
catchments_lcc.crs = lcc_crs

pft_grass = pft_wrf_reproj['GRASS-NAT']

# plot
f,ax = plt.subplots()
pft_grass.plot(ax=ax)
catchments_lcc.boundary.plot(ax=ax, color="red")
plt.show()

#%% change order of layer (first all veg PFTs, then all non-veg PFTs)
# PFTs after Harper, K. et al. (2023)

TREES_BE = np.array(pft_frac['TREES-BE']) #1 Broadleaf evergreen trees
TREES_BD = np.array(pft_frac['TREES-BD']) #2 Broadleaf deciduous trees
TREES_NE = np.array(pft_frac['TREES-NE']) #3 Needle-leaf evergreen trees
TREES_ND = np.array(pft_frac['TREES-ND']) #4 Needle-leaf deciduous trees
GRASS_NAT = np.array(pft_frac['GRASS-NAT']) #5 Grass Natural
GRASS_MAN = np.array(pft_frac['GRASS-MAN']) #6 Grass Managed (herbaceous crop land)
SHRUBS_BE = np.array(pft_frac['SHRUBS-BE']) #7 Broadleaf evergreen shrubs
SHRUBS_BD = np.array(pft_frac['SHRUBS-BD']) #8 Needle-leaf evergreen shrubs
SHRUBS_NE = np.array(pft_frac['SHRUBS-NE']) #9 Broadleaf deciduous shrubs
SHRUBS_ND = np.array(pft_frac['SHRUBS-ND']) #10 Needle-leaf deciduous shrubs
BUILT = np.array(pft_frac.BUILT) #11 Urban
WATER_INLAND = np.array(pft_frac.WATER_INLAND) #12 Inland Water
BARE = np.array(pft_frac.BARE) #13 Bare Soil
SNOWICE = np.array(pft_frac.SNOWICE) #14 Ice and Snow

# join arrays in correct order
layer_order = [TREES_BE, TREES_BD, TREES_NE, TREES_ND, GRASS_NAT, GRASS_MAN, SHRUBS_BE, SHRUBS_BD, SHRUBS_NE, SHRUBS_ND, BUILT, WATER_INLAND, BARE, SNOWICE]

# stack layers
layer_stack = np.stack(layer_order, axis=0)

#%% check if sum of all fractions adds up to 1
check_sum1 = np.sum(layer_stack, axis=0)
np.unique(check_sum1)

#%% normalize by dividing by the sum of all layer to make sure that PFTs add up to 1
# Sum along the layers axis to get the total sum for each pixel
sum_of_layers = np.sum(layer_stack, axis=0)
normalized_layer_stack = np.where(sum_of_layers != 0, layer_stack / sum_of_layers, 0)

check_sum2 = np.sum(normalized_layer_stack, axis=0)
np.unique(check_sum2) # why does it give me more than 2 values, and not only 0 and 1?

#%% set SNOWICE == BARE (soil) in pft_param.nml file

#%% create netCDF
# get lat lon values
ds_lat = pft_frac['lat'] 
ds_lon = pft_frac['lon']
# get dimensions
lat_dim = np.array(ds_lat.lat)
lon_dim = np.array(ds_lon.lon)
#get coordinates
lat = np.array(ds_lat)
lon = np.array(ds_lon)

#%% create netCDF 
Nc_img = xr.Dataset(
    coords = {'lon': lon_dim, 'lat': lat_dim, 'type': np.arange(1,len(layer_order)+1)},
    data_vars = {
    'frac': xr.DataArray(
        data = normalized_layer_stack, 
        dims = ['type','lat', 'lon']
        ),
    },
   )
Nc_img = Nc_img.transpose('type', 'lat', 'lon') # in case it need re-ordering to be displayed with ncview, change order of type, lat and lon

#%% save land cover fraction information as netCDF file to be used in JULES
Nc_img.to_netcdf('~/JULES_preprocessing/Output/netcdf/Output/netcdf/jules_pft_14_lcc_dar_d04.nc')




