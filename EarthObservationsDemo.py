# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 17:41:33 2026

@author: Refilwe Kai-Sikhakhane
"""

# Creating a subset of the data

import xarray as xr
import numpy as np

ds = xr.open_dataset(
    "S5P_RPRO_L2__NO2____20210519T120445_20210519T134615_18644_03_020400_20221107T225019.nc",
    group="PRODUCT"
)
print(ds)

print("\nVariables:")
print(list(ds.data_vars))

print("\nDimensions:")
print(ds.dims)

vars_keep = [
    "nitrogendioxide_tropospheric_column",
    "latitude",
    "longitude"
]

ds = ds[vars_keep]

# South Africa bounding box
ds = ds.where(
    (ds.latitude >= -35) & (ds.latitude <= 0) &
    (ds.longitude >= 0) & (ds.longitude <= 45),
    drop=True
)

# Save subset
ds.to_netcdf("NO2_SouthAfrica_demo.nc")

print("Subset file created: NO2_SouthAfrica_demo.nc")

# --------------------------------------------------
# Part 2
# --------------------------------------------------

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


FILE_NAME = "NO2_SouthAfrica_demo.nc"
MAKE_MAP = True

ds = xr.open_dataset(FILE_NAME)
print(ds)

# If using full path (Your data file is somewhere else on your computer)
#FILE_PATH = r"C:\Users\...\...\NO2_SouthAfrica_demo.nc" #Windows users  
#FILE_PATH = "/file/path/.../NO2_SouthAfrica_demo.nc" #For Linux or Mac users
#ds = xr.open_dataset(FILE_PATH)
#print(ds)

# --------------------------------------------------
# Select variable based on product type
# --------------------------------------------------
if "NO2" in FILE_NAME:
    var_name = "nitrogendioxide_tropospheric_column"
    units_label = "mol m$^{-2}$"
elif "AER_AI" in FILE_NAME:
    var_name = "aerosol_index_354_388"
    units_label = "Aerosol Index"
else:
    raise ValueError("Unknown product type")

data = ds[var_name]

# --------------------------------------------------
# Extract geolocation
# --------------------------------------------------
lat  = ds["latitude"]
lon  = ds["longitude"]
data = ds[var_name]

# --------------------------------------------------
# Mask fill values automatically
# --------------------------------------------------
#data = data.where(data != data.attrs["_FillValue"])
data = data.where(np.isfinite(data))

# --------------------------------------------------
# Basic statistics (ignoring missing values)
# --------------------------------------------------
print(f"Mean:   {float(data.mean()):.2e}")
print(f"StdDev: {float(data.std()):.2e}")
print(f"Median: {float(data.median()):.2e}")

print(f"Latitude range:  {float(lat.min()):.2f} to {float(lat.max()):.2f}")
print(f"Longitude range: {float(lon.min()):.2f} to {float(lon.max()):.2f}")

# --------------------------------------------------
# Plot
# --------------------------------------------------
if MAKE_MAP:
    # Create figure
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Set extent to South Africa
    ax.set_extent([0, 45, -35, 0], crs=ccrs.PlateCarree())

    # Clearer geographic context
    ax.add_feature(cfeature.COASTLINE, linewidth=1.2, edgecolor="black")
    ax.add_feature(cfeature.BORDERS, linewidth=1.2, edgecolor="black")

    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.4,
        color="gray",
        alpha=0.6,
        linestyle="--"
    )
    gl.top_labels = False
    gl.right_labels = False

    # Mask missing values and clip negatives
    data_clean = data.where(data.notnull()).clip(min=0)

    # Determine color scale
    vmax = float(data_clean.max()) * 0.05  # adjust scaling if needed

    try:
        # Attempt pcolormesh (gridded case)
        mesh = ax.pcolormesh(
            lon, lat, data_clean,
            cmap="viridis",
            vmin=0,
            vmax=vmax,
            transform=ccrs.PlateCarree()
        )
        cb = plt.colorbar(mesh, ax=ax, orientation="vertical", pad=0.02)
        cb.set_label(units_label)
        ax.set_title("Sentinel-5P NO₂ Tropospheric Column")

    except ValueError:
        # Fallback: swath scatter plot
        sc = ax.scatter(
            lon, lat, c=data_clean,
            s=1,
            cmap="viridis",
            vmin=0,
            vmax=vmax,
            transform=ccrs.PlateCarree()
        )
        cb = plt.colorbar(sc, ax=ax, orientation="vertical", pad=0.02)
        cb.set_label(units_label)
        ax.set_title("Sentinel-5P NO₂ Swath data")

    plt.show()

