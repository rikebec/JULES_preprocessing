import xarray as xr
import numpy as np
import os
import xesmf as xe
import math
import pandas as pd

# --- User-defined variables ---
catchment = 'Vilcanota'
# Note: Adjust these paths to your actual file locations
COSIPY_dir = '/home/rike/JULES_postprocessing/COSIPY_output/'
JULES_dir = '/home/rike/dar_data/d04/'+catchment+'/ancil/'
cosipy_file = os.path.join(COSIPY_dir, "WRF_ERA_Vilca_500m_runoff_20000101-20001231.nc")
jules_file = os.path.join(JULES_dir, "LandFrac_"+catchment+".nc")

# --- Data Loading ---
# Load high-resolution glacier data
glacier_ds = xr.open_dataset(cosipy_file)
# Load low-resolution grid data for target extent
grid_ds = xr.open_dataset(jules_file)

# --- Data Processing ---
# Combine melt and precipitation to get total runoff
glacier_runoff = glacier_ds["Q"] + glacier_ds["RAIN"]
glacier_runoff.name = "glacier_runoff"
print('Melt and rain combined.')

# --- Grid Extension ---
# Define the target grid using the extent of the coarse grid
# and the resolution of the fine grid.

# 1. Get the target extent from the low-resolution grid
lat_min, lat_max = grid_ds['lat'].min().item(), grid_ds['lat'].max().item()
lon_min, lon_max = grid_ds['lon'].min().item(), grid_ds['lon'].max().item()

# 2. Get the resolution from the high-resolution glacier data
# Ensure coordinates are sorted to correctly calculate the difference.
glacier_runoff = glacier_runoff.sortby('lat').sortby('lon')
lat_res = np.abs(np.diff(glacier_runoff['lat'])[0])
lon_res = np.abs(np.diff(glacier_runoff['lon'])[0])

# 3. Create the new, extended high-resolution coordinate arrays
# Using np.arange is precise for creating sequences with a fixed step.
new_lats = np.arange(lat_min, lat_max + lat_res/2, lat_res)
new_lons = np.arange(lon_min, lon_max + lon_res/2, lon_res)

# 4. Use reindex to place the original data onto the new extended grid
# This is the key step. It aligns data based on coordinate labels and is the
# idiomatic way to perform this operation in xarray.
# It fills areas outside the original data's extent with NaN.
print("Extending grid using reindex...")
glacier_runoff_extended = glacier_runoff.reindex(
    {"lat": new_lats, "lon": new_lons},
    method="nearest",
    tolerance=lat_res # Use a tolerance to handle floating point inaccuracies
)

# --- Verification and Output ---
# Verify that no duplicate coordinates were created
lat_dupes = np.any(xr.DataArray(new_lats).to_series().duplicated())
lon_dupes = np.any(xr.DataArray(new_lons).to_series().duplicated())
print(f"Duplicate latitudes found in new grid: {lat_dupes}")
print(f"Duplicate longitudes found in new grid: {lon_dupes}")

if not lat_dupes and not lon_dupes:
    # Write the final extended data to a new NetCDF file
    output_filename = "Extended_glacier_runoff.nc"
    glacier_runoff_extended.to_netcdf(output_filename)
    print(f"Successfully created '{output_filename}' with the new extent.")
else:
    print("Duplicates found. Output file not written.")


# --- get fraction of glacier coverage --- #

print("Calculating fraction of glacier coverage...")

# 1. Create a binary mask: 1 where glaciers exist, 0 otherwise.
# runoff grid is considered a glacier pixel.
binary_mask = xr.where(glacier_runoff_extended.notnull(), 1, 0)
binary_mask.name = "glacier_mask"

print("Binary mask created.")
'''
# load for debugging and to skip first steps if already run.
glacier_runoff_extended = xr.open_dataset('Extended_glacier_runoff.nc')
glacier_runoff_extended = glacier_runoff_extended['glacier_runoff']
binary_mask = xr.open_dataset('binary_mask.nc')
binary_mask = binary_mask['glacier_runoff']
grid_ds = xr.open_dataset(jules_file)
grid_ds = grid_ds["land_frac"]
'''

# --- Diagnostic Print: Check Binary Mask ---
num_glacier_pixels = int(binary_mask.sum().item())
print(f"Binary mask created. Total number of glacier pixels (value=1): {num_glacier_pixels}")
if num_glacier_pixels == 0:
    print("Warning: The binary mask is all zeros. The final output will also be all zeros.")

# 2. Define the aggregation bins from the target grid's coordinate centers.
# We calculate the edges of the coarse grid cells to use them as bins.
lat_coords = grid_ds.lat.values
lon_coords = grid_ds.lon.values
lat_diff = np.diff(lat_coords) / 2
lon_diff = np.diff(lon_coords) / 2
lat_bins = np.concatenate([[lat_coords[0] - lat_diff[0]], lat_coords[:-1] + lat_diff, [lat_coords[-1] + lat_diff[-1]]])
lon_bins = np.concatenate([[lon_coords[0] - lon_diff[0]], lon_coords[:-1] + lon_diff, [lon_coords[-1] + lon_diff[-1]]])

# 3. Group the binary mask by the bins and sum the values.
# This counts the number of high-res '1' pixels within each low-res bin.
# The `sum(skipna=False)` ensures that if all pixels in a bin are NaN, the result is 0, not NaN.
glacier_pixel_count = binary_mask.groupby_bins('lat', bins=lat_bins).sum(skipna=False).groupby_bins('lon', bins=lon_bins).sum(skipna=False)

# 4. Clean up coordinates for the final output.
# The coordinates are now intervals; we'll rename them and assign the original target grid centers.
glacier_pixel_count = glacier_pixel_count.rename({'lat_bins': 'lat', 'lon_bins': 'lon'})
glacier_pixel_count['lat'] = grid_ds.lat
glacier_pixel_count['lon'] = grid_ds.lon
glacier_pixel_count.name = "glacier_pixel_count"

# 5. Save the result.
output_count_filename = "Glacier_pixel_count_per_JULES_pixel.nc"
glacier_pixel_count.to_netcdf(output_count_filename)

print(f"\nSuccessfully created '{output_count_filename}'.")
print(f"Max pixel count in a JULES cell: {glacier_pixel_count.max().item()}")

temp = glacier_runoff_extended.copy()
temp = temp.sortby('lat').sortby('lon')
lat_res = np.abs(np.diff(temp['lat'])[0])
lon_res = np.abs(np.diff(temp['lon'])[0])
del temp

target_lat_res = np.abs(np.diff(grid_ds.lat.values)).mean()
target_lon_res = np.abs(np.diff(grid_ds.lon.values)).mean()
max_possible_pixels = (target_lat_res / lat_res) * (target_lon_res / lon_res)
print(f"\nMaximum possible high-res pixels per low-res cell (100% coverage): {max_possible_pixels:.2f}")

glacier_frac = glacier_pixel_count/math.ceil(max_possible_pixels)

glacier_frac.to_netcdf("Glacier_coverage_frac_per_JULES_pixel.nc")

# --- get total runoff from COSIPY-pixels within each JULES-pixel ---
# 3. Group the binary mask by the bins and sum the values.
# This counts the number of high-res '1' pixels within each low-res bin.
# The `sum(skipna=False)` ensures that if all pixels in a bin are NaN, the result is 0, not NaN.
glacier_runoff_sum = glacier_runoff_extended.groupby_bins('lat', bins=lat_bins).sum(skipna=True).groupby_bins('lon', bins=lon_bins).sum(skipna=True)

# 4. Clean up coordinates for the final output.
# The coordinates are now intervals; we'll rename them and assign the original target grid centers.
glacier_runoff_sum = glacier_runoff_sum.rename({'lat_bins': 'lat', 'lon_bins': 'lon'})
glacier_runoff_sum['lat'] = grid_ds.lat
glacier_runoff_sum['lon'] = grid_ds.lon
glacier_runoff_sum.name = "glacier_runoff"

#print("get total runoff from COSIPY-pixels within each JULES-pixel")
#sum_COSIPY_runoff = regridder(glacier_runoff_extended)
#sum_COSIPY_runoff.name = 'runoff_COSIPY'
glacier_runoff_sum.to_netcdf("COSIPY_runoff_per_JULES_pixel.nc")

# --- load the surface runoff from JULES output ---
jules_surface_runoff = xr.open_dataset("../output_with_wrf_ice_2010_2019/jules-vn6.1.S2.daily_hydrology.2012.2D.nc")
jules_surface_runoff = jules_surface_runoff["surf_roff"]

# add the time to glacier_frac and glacier_runoff
datetime_cosipy = pd.date_range("2012-01-01 00:00:00", periods = 8784, freq = "H")
glacier_runoff_sum["time"] = datetime_cosipy
glacier_frac["time"] = datetime_cosipy

# resample to daily
glacier_runoff_sum = glacier_runoff_sum.resample(time = "1D").sum()
glacier_frac = glacier_frac.resample(time = "1D").mean()

jules_final = ((1 - glacier_frac) * jules_surface_runoff) + glacier_runoff_sum

jules_final.to_netcdf("jules+cosipy_2012_fake.nc")
