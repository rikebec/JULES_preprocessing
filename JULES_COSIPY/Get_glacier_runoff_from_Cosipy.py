import xarray as xr
import numpy as np
import os
import xesmf as xe

# --- User-defined variables ---
catchment = 'Vilcanota'
# Note: Adjust these paths to your actual file locations
COSIPY_dir = '/home/rike/JULES_postprocessing/COSIPY_output/'
JULES_dir = '/home/rike/dar_data/d04/'+catchment+'/ancil/'
# Load low-resolution grid data for target extent
grid_ds = xr.open_dataset(jules_file)

# --- Data Processing ---
# Combine melt and precipitation to get total runoff

# --- Grid Extension ---
# Define the target grid using the extent of the coarse grid
# and the resolution of the fine grid.
# Note: Adjust these paths to your actual file locations
COSIPY_dir = '/home/rike/JULES_postprocessing/COSIPY_output/'
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
# Define the target grid using the extent of the coarse grid
# and the resolution of the fine grid.
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
print("getting fraction of glacier coverage ...")
# Set all NANs to zero and all glacier pixels to 1
binary_mask = xr.where(np.isnan(glacier_runoff_extended), 0, 1)
binary_mask.to_netcdf("binary_mask.nc")
print("binary mask created")

# resample from COSIPY-grid to JULES-grid
source_grid = xr.Dataset({
    'lat': (['lat'], glacier_runoff_extended.lat.values),
    'lon': (['lon'], glacier_runoff_extended.lon.values)
})

target_grid = xr.Dataset({
    'lat': (['lat'], grid_ds.lat.values),
    'lon': (['lon'], grid_ds.lon.values)
})

regridder = xe.Regridder(source_grid, target_grid, method='conservative', reuse_weights=False)
aggregated = regridder(binary_mask)
# get fraction of COSIPY-pixel coverage within each JULES-pixel
normalized = aggregated / aggregated.max()
normalized.to_netcdf("Glacier_coverage_frac_per_JULES_pixel.nc")

# --- get total runoff from COSIPY-pixels within each JULES-pixel --- #
print("get total runoff from COSIPY-pixels within each JULES-pixel")
sum_COSIPY_runoff = regridder(glacier_runoff_extended)
sum_COSIPY_runoff.to_netcdf("COSIPY_runoff_per_JULES_pixel.nc")
