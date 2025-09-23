
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

STD_LEVELS = [800, 850]

STD_HOURS = []

lat = 35.512219
lon = -108.822015

warnings.filterwarnings("ignore")

save_dir = "/Users/ryanpurciel/Documents/case/"

path = "/Users/ryanpurciel/Documents/case/Data/hrrr/"

files = sorted(glob.glob(path + "*.grib2"))

print("Starting...")

meteogram_df = pd.DataFrame(columns = ("time_UTC", "pressure", "temp_c", "dpt_c", "rh_pct", "ws_kt", "wd_deg", "wg_kt", "precip_type", "precip_rate_kg_m2s", "frozen_precip_pct"))

for file in files:

	print("Selected File ", file)

	sfc_inst_data = xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'surface', 'stepType': 'instant'})
	m2_data = xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 2})
	m10_data = xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 10})
	
	# Becuase the grib reader reads longitude from 0-360 and not -180-180
	# we have to adjust the `lon`.
	if lon < 0:
		lon += 360
	
	# Find the nearest neightbor from the maximum absolute differences...
	abslat = np.abs(sfc_inst_data.latitude-lat)
	abslon = np.abs(sfc_inst_data.longitude-lon)

	c = np.maximum(abslon, abslat)
	
	# Get the index of the nearest point
	# (Again, the values are backwards you might expect because the 
	# coordinates are in (y, x) order.)
	([idx_y], [idx_x]) = np.where(c == np.min(c))
	
	#Now that we have the index location of the nearest point to the requested lat/lon, 
	#we can use the select function sel to get all the data in the height coordinate at a single point in the dataset.
	
	point_sfc_inst = sfc_inst_data.sel(y=idx_y, x=idx_x)
	point_2m = m2_data.sel(y=idx_y, x=idx_x)
	point_10m = m10_data.sel(y=idx_y, x=idx_x)

	print(f"SFC pt: {point_sfc_inst.latitude.values}, {point_sfc_inst.longitude.values-360}")
	print(f"2m pt: {point_2m.latitude.values}, {point_2m.longitude.values-360}")
	print(f"10m pt: {point_10m.latitude.values}, {point_10m.longitude.values-360}")

	time_utc = pd.to_datetime(point_sfc_inst.time.values)
	print(time_utc)
	
	pressure_pa = units.Quantity(point_sfc_inst.sp.values, 'pascal')
	pressure_hPa = pressure_pa.to(units.hectopascal)
	print(pressure_hPa)
	
	temp_k = units.Quantity(point_2m.pt.values, 'K')
	temp_c = temp_k.to(units.celsius)
	print(temp_c)

	dpt_k = units.Quantity(point_2m.d2m.values, 'K')
	dpt_c = dpt_k.to(units.celsius)
	print(dpt_c)

	rh = units.Quantity(point_2m.r2.values, '%')
	print(rh)

	u = units.Quantity(point_10m.u10.values, 'm/s')
	v = units.Quantity(point_10m.u10.values, 'm/s')
	ws_ms = mpcalc.wind_speed(u, v)
	ws_kts = ws_ms.to(units.knots)

	wd_deg = mpcalc.wind_direction(u, v)

	print(ws_kts)
	print(wd_deg)

	wg_ms = units.Quantity(point_sfc_inst.gust.values, 'm/s')
	wg_kts = wg_ms.to(units.knots)
	print(wg_kts)

	precip_rate = units.Quantity(point_sfc_inst.prate.values, "kg m**-2 s**-1")
	print(precip_rate)

	pct_frz_precip = units.Quantity(point_sfc_inst.cpofp.values, "%")
	print(pct_frz_precip)

	cat_ra = point_sfc_inst.crain.values
	cat_fzrn = point_sfc_inst.cfrzr.values
	cat_ip = point_sfc_inst.cicep.values
	cat_sn = point_sfc_inst.csnow.values

	model_ptypes_name = ["Rain", "Freezing Rain", "Ice Pellets", "Snow"]
	model_ptypes_present = [cat_ra, cat_fzrn, cat_ip, cat_sn]

	precip_type = "None"
	for pval, pname in zip(model_ptypes_present, model_ptypes_name):
		print(f"{pname}: {pval}")
		if pval == 1:
			precip_type += f"{pname}, "

	if precip_type != "None":
		incl_precip_type = precip_type[:-1].replace("None", "")
		incl_pct_frz_precip = pct_frz_precip.magnitude
		incl_precip_rate = precip_rate.magnitude
	else:
		incl_precip_type = precip_type
		incl_pct_frz_precip = 0
		incl_precip_rate = 0

	print(incl_precip_type, incl_precip_rate, incl_pct_frz_precip)
	
	meteogram_df.loc[len(meteogram_df)] = [time_utc.strftime("%Y-%m-%d %H:%M:%S"),
										   pressure_hPa.magnitude,
										   temp_c.magnitude,
										   dpt_c.magnitude,
										   rh.magnitude,
										   ws_kts.magnitude,
										   wd_deg.magnitude,
										   wg_kts.magnitude,
										   incl_precip_type,
										   incl_precip_rate,
										   incl_pct_frz_precip,
										   ]
	print(meteogram_df)

meteogram_df.to_csv(os.path.join(save_dir, "meteogram_data.csv"), index=False)

print("DONE!!!!!!!!!!!!!!!!! FINALLY!!!!!!!!!!!!!!")



