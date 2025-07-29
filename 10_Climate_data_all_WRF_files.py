#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 10:13:37 2024

@author: rikebecker

This script pre-processes the WRF climate data, created for the Deplete and Retreat project, 
to make the files readable by the JULES model.

The temporal resolution of the WRF data is hourly!

!!! still work in progress. To be run later on JASMIN server to process data stored in the Deplete and Retreat group workspace.

"""

#%% load libraries
import xarray as xr
import os
#import matplotlib.pyplot as plt
#import rioxarray as rio
#import rasterio
import numpy as np
import salem
#import pandas as pd
#import geopandas as gpd

#%% set the path to the directory where the input and output data is stored
I_directory = '.../Input/climate/d03'
O_directory = '.../Output/climate/d03'

#%% load WRF output files
file_list = os.listdir(I_directory) #list all files and put into loop


for file_name in file_list:
    img_name = I_directory+'/'+file_name
    date = file_name[-20:]

    # read the image 
    ds = salem.open_wrf_dataset(img_name, decode_times=False)

    ## ===================================================================== ##
    ##                          Temperature (T2) in K                        ##
    ## ===================================================================== ##
    temp = ds.T2 
    temp = temp.drop_vars(["lat", "lon", "xtime"])
    temp = temp.rename({"south_north":"lat","west_east":"lon"})
    temp.name = 't'

    #add attributes
    temp.attrs['units'] = 'K'
    temp.attrs['standard_name'] = 't'
    temp.attrs['long_name'] = 'Air temperature at 2m'

    #write out and save netcdf in output directory
    temp.to_netcdf(O_directory+"/t"+date+".nc")
 
    ## ===================================================================== ##
    ##                 Precipitation (RAINC + RAINNC) in mm                  ##
    ## ===================================================================== ##

    rainnc = ds.RAINNC
    rainc = ds.RAINC #likely to be zero
    rain_cumulative = rainc+rainnc

    #convert from cumulative to non-cumulative data
    rain_hourly = rain_cumulative.diff(dim='time')
    initial_values = xr.zeros_like(rain_cumulative.isel(time=0))  # to keep original number of time steps
    rain_hourly = xr.concat([initial_values, rain_hourly], dim='time')
    # Ensure the dimension order remains the same
    rain_hourly = rain_hourly.assign_coords(time=rain_cumulative['time'])
    rain_hourly = rain_hourly.transpose(*rain_cumulative.dims)

    #convert from mm/hour to kg m-2 s-1
    #conversion_factor = 1/1800 #for precip data in mm/30min
    conversion_factor = 1/3600 #for precip data in mm/hour
    
    precip = rain_hourly*conversion_factor
    precip = precip.drop_vars(["lat", "lon", "xtime"])
    precip = precip.rename({"south_north":"lat","west_east":"lon"})

    #add attributes
    precip.attrs['units'] = 'kg m-2 s-1'
    precip.attrs['standard_name'] = 'precip'
    precip.attrs['long_name'] = 'hourly total precipitation'
    #write out and save netcdf in output directory
    precip.to_netcdf(O_directory+"/precip"+date+".nc")

    ## ===================================================================== ##
    ##                     Specific Humidity (Q2) in kg kg-1                 ##
    ## ===================================================================== ##

    hum_spec = ds.Q2 
    hum_spec = hum_spec.drop_vars(["lat", "lon", "xtime"])
    hum_spec = hum_spec.rename({"south_north":"lat","west_east":"lon"})

    #add attributes
    hum_spec.attrs['units'] = 'kg kg-1'
    hum_spec.attrs['standard_name'] = 'q'
    hum_spec.attrs['long_name'] = 'specific humidity at 2m'
    #write out and save netcdf in output directory
    hum_spec.to_netcdf(O_directory+"/q"+date+".nc")

    ## ===================================================================== ##
    ##                    Wind speed (U10 and V10) in m s-1                  ##
    ## ===================================================================== ##

    wind_u = ds.U10
    wind_v = ds.V10
    # get total wind speed (note! JULES can also handle U and V variables!)
    wind = np.sqrt(wind_u**2 + wind_v**2)
    
    wind = wind.drop_vars(["lat", "lon", "xtime"])
    wind = wind.rename({"south_north":"lat","west_east":"lon"})

    #add attributes
    wind.attrs['units'] = 'm s-1'
    wind.attrs['standard_name'] = 'wind'
    wind.attrs['long_name'] = 'wind speed at 2m'
    #write out and save netcdf in output directory
    wind.to_netcdf(O_directory+"/wind"+date+".nc")

    ## ===================================================================== ##
    ##          Shortwave downwelling radiation (SWDOWN) in W m-2           ##
    ## ===================================================================== ##
    
    sw_down = ds.SWDOWN
    sw_down = sw_down.drop_vars(["lat", "lon", "xtime"])
    sw_down = sw_down.rename({"south_north":"lat","west_east":"lon"})

    #add attributes
    sw_down.attrs['units'] = 'W m-2'
    sw_down.attrs['standard_name'] = 'sw_down'
    sw_down.attrs['long_name'] = 'Downward short wave flux at ground surface'
    #write out and save netcdf in output directory
    sw_down.to_netcdf(O_directory+"/sw_down"+date+".nc")

 
    ## ===================================================================== ##
    ##              Longwave downwelling radiation (GLW) in W m-2           ##
    ## ===================================================================== ##

    lw_down = ds.GLW
    lw_down = lw_down.drop_vars(["lat", "lon", "xtime"])
    lw_down = lw_down.rename({"south_north":"lat","west_east":"lon"})

    #add attributes
    lw_down.attrs['units'] = 'W m-2'
    lw_down.attrs['standard_name'] = 'lw_down'
    lw_down.attrs['long_name'] = 'Downward long wave flux at ground surface'
    #write out and save netcdf in output directory
    lw_down.to_netcdf(O_directory+"/lw_down"+date+".nc")

    ## ===================================================================== ##
    ##                     Surface pressure (PSFC) in Pa                     ##
    ## ===================================================================== ##

    pstar = ds.PSFC
    pstar = pstar.drop_vars(["lat", "lon", "xtime"])
    pstar = pstar.rename({"south_north":"lat","west_east":"lon"})

    #add attributes
    pstar.attrs['units'] = 'Pa'
    pstar.attrs['standard_name'] = 'pstar'
    pstar.attrs['long_name'] = 'Pressure at surface'
    #write out and save netcdf in output directory
    pstar.to_netcdf(O_directory+"/pstar"+date+".nc")


#%% join all files for each variable
# consider using CDO for a faster joining of the single files

list1 = ["/t_*", "/precip_*", "/q_*", "/wind_*", "/sw_*", "/lw_*", "/pstar_*"]
list2 = ["/t.nc", "/precip.nc", "/q.nc", "/wind.nc", "/swdown.nc", "/lwdown.nc", "/pstar.nc"]

for var1, var2 in zip(list1, list2):
    print("processing: "+var2)   
    ds = xr.open_mfdataset(O_directory+var1)
    ds.compute()
    ds.to_netcdf(O_directory+var2)
    # Delete single files
    files = glob.glob(os.path.join(O_directory+var1))
    for file in files:
        os.remove(file)
