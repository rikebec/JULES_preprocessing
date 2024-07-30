<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 09:47:57 2024

@author: Rike Becker

Description: cut shapefile into Deplete and Retreat domains 03 
after re-projecting it to Lambert Conformal Conic (LCC) projection
to match the WRF-projection
"""
#%% load libraries
import geopandas as gpd
from shapely.geometry import box
import matplotlib.pyplot as plt
import salem # for the transformation of the coordinate systems
import numpy as np

#%% read the shapefile and get extent
fp = "~/JULES_preprocessing/Input/shapefiles/Basins_selected.shp"
sf = gpd.read_file(fp, crs="epsg:4326") # read shapefile
sf.plot()
sf.crs

#%% read in the land_frac file to get its coordinate system for later transformation
wrf_land_frac = "~/JULES_preprocessing/Output/netcdf/jules_land_frac_wrf_dar_d03.nc" 
ds = salem.open_wrf_dataset(wrf_land_frac)
ds.salem.grid # check crs

#%% get extent of DAR domain 03 (catchments in Chile and Argentina)
lon_max = np.array(ds.longitude).max()
lon_min = np.array(ds.longitude).min()
lat_max = np.array(ds.latitude).max()
lat_min = np.array(ds.latitude).min()

# Create a polygon for domain
polygon = box(lon_min, lat_min, lon_max, lat_max)
domain = gpd.GeoDataFrame([1], geometry=[polygon], crs=sf.crs)

#%% plot to check extent of domains
fig, (ax1) = plt.subplots(figsize=(12, 8))
sf.plot(ax=ax1)
domain.boundary.plot(ax=ax1, color="red")
plt.show()

#%% >>> clip the initial shapefile using domain_03 
catchments_clipped = sf.clip(domain)

#%% transform the shapefile to Lambert Conformal Conic (WRF-projection)
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

catchments_clipped_proj = catchments_clipped.to_crs(lcc_crs)
catchments_clipped_proj.crs

#%% plot to check extent of domains
fig,ax=plt.subplots()
ds.land_frac.plot(ax=ax) 
catchments_clipped_proj.plot(ax=ax, facecolor='none', edgecolor = 'k') 

#%% save as new shapefiles
catchments_clipped_proj.to_file("~/JULES_preprocessing/Output/shapefiles/Domain03_catchments_lcc.shp")
=======
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 09:47:57 2024

@author: Rike Becker

Description: cut shapefile into Deplete and Retreat domains 03 
after re-projecting it to Lambert Conformal Conic (LCC) projection
to match the WRF-projection
"""
#%% load libraries
import geopandas as gpd
from shapely.geometry import box
import matplotlib.pyplot as plt
import salem # for the transformation of the coordinate systems
import numpy as np

#%% read the shapefile and get extent
fp = "~/JULES_preprocessing/Input/shapefiles/Basins_selected.shp"
sf = gpd.read_file(fp, crs="epsg:4326") # read shapefile
sf.plot()
sf.crs

#%% read in the land_frac file to get its coordinate system for later transformation
wrf_land_frac = "~/JULES_preprocessing/Output/netcdf/jules_land_frac_wrf_dar_d03.nc" 
ds = salem.open_wrf_dataset(wrf_land_frac)
ds.salem.grid # check crs

#%% get extent of DAR domain 03 (catchments in Chile and Argentina)
lon_max = np.array(ds.longitude).max()
lon_min = np.array(ds.longitude).min()
lat_max = np.array(ds.latitude).max()
lat_min = np.array(ds.latitude).min()

# Create a polygon for domain
polygon = box(lon_min, lat_min, lon_max, lat_max)
domain = gpd.GeoDataFrame([1], geometry=[polygon], crs=sf.crs)

#%% plot to check extend of domains
fig, (ax1) = plt.subplots(figsize=(12, 8))
sf.plot(ax=ax1)
domain.boundary.plot(ax=ax1, color="red")
plt.show()

#%% >>> clip the initial shapefile using domain_03 
catchments_clipped = sf.clip(domain)

#%% transform the shapefile to Lambert Conformal Conic (WRF-projection)
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

catchments_clipped_proj = catchments_clipped.to_crs(lcc_crs)
catchments_clipped_proj.crs

#%% plot to check extend of domains
fig,ax=plt.subplots()
ds.land_frac.plot(ax=ax) 
catchments_clipped_proj.plot(ax=ax, facecolor='none', edgecolor = 'k') 

#%% save as new shapefiles
catchments_clipped_proj.to_file("~/JULES_preprocessing/Output/shapefiles/Domain03_catchments_lcc.shp")
>>>>>>> 4bcd258 (new and updated scripts)
