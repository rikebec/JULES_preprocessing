# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 16:46:20 2024

@author: Rike Becker

Description: script to create the land fraction input for JULES and mask out
the DAR catchments
"""
#%% load required libraries
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import mapping

#%% set the paths
# to land fraction file 
land_frac_img = "~/JULES_preprocessing/Output/netcdf/jules_land_frac_lcc_dar_d04.nc" 
# to the shapefile
catchment_file = "~/JULES_preprocessing/Output/shapefiles/Domain04_catchments_lcc.shp"

#%% read the images
land_frac_nc = xr.open_dataset(land_frac_img) # open netCDF
land_frac = land_frac_nc.land_frac
catchments = gpd.read_file(catchment_file) # read shapefile

#%% assign crs (lcc) to catchment shapefile
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
catchments.crs = lcc_crs

#%% assign crs to netCDF land fraction map
land_frac.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
land_frac.rio.write_crs(lcc_crs, inplace=True)

#%% clip land_frac with catchment shapefile (domain)
land_frac_clipped = land_frac.rio.clip(catchments.geometry.apply(mapping), catchments.crs, all_touched=True, drop=False)

#%% set all NAN value to zero
land_frac_clipped = land_frac_clipped.fillna(0)

#%% plot
f,ax = plt.subplots()
land_frac_clipped.plot(ax=ax)
#catchments.boundary.plot(ax=ax, color="k")
plt.show()

#%% Replace old land mask with new land mask
land_frac_nc['land_frac']=land_frac_clipped

#%% save as clipped land frac
land_frac_nc.to_netcdf("~/JULES_preprocessing/Output/netcdf/jules_land_frac_lcc_dar_d04_clipped.nc")


