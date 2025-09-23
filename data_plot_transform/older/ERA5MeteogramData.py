
import numpy as np
import xarray as xr
from metpy.units import units
import metpy.interpolate as metpyinterp
from metpy.interpolate import cross_section
import metpy.calc as mpcalc
import metpy
import pandas as pd
from pandas import read_csv
import cfgrib
import sys, os
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors, to_rgba
import matplotlib.colors as mcolors
import cartopy.crs as crs
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
import cartopy
import matplotlib.ticker as mticker
import netCDF4
import math
import warnings
import glob
import statistics

lat = 3.501691
lon = 36.040758

warnings.filterwarnings("ignore")

save_dir = "/Users/rpurciel/Documents/case/2m Air Temp/"

path = "/Users/rpurciel/Documents/Bcase/ERA5/data/"

files = sorted(glob.glob(path + "*SFC_ERA5.grib"))

print("Starting...")

meteogram_df = pd.DataFrame(columns = ("time_UTC", "pressure_inhg", "temp_f", "dpt_f", "ws_kts", "wd_deg", "wg_kts"))

for file in files:

	print("Selected File ", file)

	sfc_data = xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'surface'})
	
	point_sfc = sfc_data.sel(longitude=lon, latitude=lat, method= 'nearest')

	extent = [lon-5, lon+5, lat-5, lat+5]

	fig = plt.figure(figsize=(22,16))
	ax = plt.axes(projection = ccrs.PlateCarree())
	ax.set_extent(extent, crs=ccrs.PlateCarree())
	states = NaturalEarthFeature(category="cultural", scale="50m",
	                                      facecolor="none",
	                                      name="admin_1_states_provinces")
	ax.add_feature(states, linewidth=1.0, edgecolor="black")
	ax.coastlines('50m', linewidth=1.5)
	ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), linestyle='-', linewidth=0.5, alpha=1,edgecolor='blue',facecolor='none')
	ax.add_feature(cfeature.BORDERS, linewidth=1.5)

	ax.set_xlabel('Latitude')
	ax.set_ylabel('Longitude')

	gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
	gl.xlabels_top = False
	gl.ylabels_left = True
	gl.ylabels_right= False
	gl.xlines = True
	gl.ylines = True
	gl.xlocator = mticker.FixedLocator(np.arange(extent[0], extent[1], 2))
	gl.ylocator = mticker.FixedLocator(np.arange(extent[3], extent[2], 1))

	plt.plot(lon, lat, marker='X', markersize=8, markerfacecolor='greenyellow', markeredgecolor='black', transform=ccrs.PlateCarree(), zorder=15)

	time_utc = pd.to_datetime(point_sfc.time.values)
	print(time_utc)

	plt.savefig(os.path.join(save_dir, f"verification_{time_utc}.png"))
	
	pressure_pa = units.Quantity(point_sfc.sp.values, 'pascal')
	pressure_hPa = pressure_pa.to(units.hectopascal)
	pressure_inHg = pressure_pa.to(units.inHg)
	print(pressure_hPa, pressure_inHg)
	
	temp_k = units.Quantity(point_sfc.t2m.values, 'K')
	temp_c = temp_k.to(units.celsius)
	temp_f = temp_k.to(units.fahrenheit)
	print(temp_c, temp_f)

	dpt_k = units.Quantity(point_sfc.d2m.values, 'K')
	dpt_c = dpt_k.to(units.celsius)
	dpt_f = dpt_k.to(units.fahrenheit)
	print(dpt_c, dpt_f)

	u = units.Quantity(point_sfc.u10.values, 'm/s')
	v = units.Quantity(point_sfc.u10.values, 'm/s')
	ws_ms = mpcalc.wind_speed(u, v)
	ws_kts = ws_ms.to(units.knots)

	wd_deg = mpcalc.wind_direction(u, v)

	print(ws_kts)
	print(wd_deg)
	
	meteogram_df.loc[len(meteogram_df)] = [time_utc.strftime("%Y-%m-%d %H:%M:%S"),
										   pressure_inHg.magnitude,
										   temp_f.magnitude,
										   dpt_f.magnitude,
										   ws_kts.magnitude,
										   wd_deg.magnitude,
										   -999
										   ]
	print(meteogram_df)

meteogram_df.to_csv(os.path.join(save_dir, "meteogram_data.csv"), index=False)

print("DONE!!!!!!!!!!!!!!!!! FINALLY!!!!!!!!!!!!!!")



