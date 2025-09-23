
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
import matplotlib.patheffects as PathEffects
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

from __internal_funcs import plot_towns, draw_logo

STD_LEVELS = [9999]

PLOT_SFC_WIND_GUSTS = True

STD_HOURS = []

PRODUCT = "Wind"

idx = 8

wind_display = "Barbs"

warnings.filterwarnings("ignore")

save_dir = "/Users/rpurciel/Documents/case/HRRR/Sfc Gusts/"

path = "/Users/rpurciel/Documents/case/HRRR/Data/"

files = sorted(glob.glob(path + "*.grib2"))

print("Starting...")

for file in files:

	print("Selected File ", file)

	raw_data = []
	raw_data += [xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'isobaricInhPa'})]
	raw_data = xr.merge(raw_data)
	sfc_data = xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'surface', 'stepType': 'instant'})
	m2_data = xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 2})
	m10_data = xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 10})
	

	num_times = 1
	num_levels = len(raw_data.isobaricInhPa)

	for time in range(0, num_times, 1):

		timesel_data = raw_data

		sel_time = str(timesel_data['time'].values)

		#print(timesel_data['isobaricInhPa'].values)

		print("		Time: ", sel_time)

		for level in STD_LEVELS:

			if level == 9999:
				sel_level = 'SFC'
				data = sfc_data
			else:
				sel_level = str(level)
				data = timesel_data.sel(isobaricInhPa=level)

			print("			Level: ", sel_level)

			# extent = [250,359,0,65]

			# nd1 = levelsel_data.where(levelsel_data.latitude <= 65, drop=True) 
			# nd2 = nd1.where(levelsel_data.latitude >= -0, drop=True)
			# nd3 = nd2.where(levelsel_data.longitude <= 359, drop=True)
			# nd4 = nd3.where(levelsel_data.longitude >= 250, drop=True)

			# data = nd4

			lat = data.variables['latitude'][:]
			lon = data.variables['longitude'][:]

			if level == 9999: #sfc flag
				u = m10_data.variables['u10'][:]
				v = m10_data.variables['v10'][:]
				ukt = units.Quantity(u.values, 'm/s').to(units.knots)
				vkt = units.Quantity(v.values, 'm/s').to(units.knots)

				t = m2_data.variables['t2m'][:]
				t_c = units.Quantity(t.values, 'K').to(units.celsius)

				z = sfc_data.variables['orog'][:]
				gpm = units.Quantity(z.values, 'm')

				if PLOT_SFC_WIND_GUSTS:
					ws_ms = sfc_data.variables['gust'][:]
					ws = units.Quantity(ws_ms.values, 'm/s').to(units.knots)
					title_desc = "10m Wind Speed (kts, barbs)\nand Instantaneous Wind Gust (kts, shaded)" #, 2m Temperature (°C, colored contours),  and Terrain Height (m, contours)
				else:
					ws = mpcalc.wind_speed(m10_data['u'], m10_data['v']).metpy.convert_units('knots')
					title_desc = "10m Wind Speed (kts, shaded), 2m Temperature (°C, colored contours),\nand Terrain Height (m, contours)"

				gpm_center = round((gpm.magnitude.min() + gpm.magnitude.max())/2, 10)
				gpm_levels = np.arange(gpm.magnitude.min(), gpm.magnitude.max(), 250)

				t_center = round((t_c.magnitude.min() + t_c.magnitude.max())/2, 10)
				t_levels = np.arange(t_center - 20, t_center + 20, 8)

				title_info = f'Surface Level\nValid at {np.datetime_as_string(data["valid_time"].values, timezone="UTC")[:-11].replace("T", " ")}'

			else:
				u = data.variables['u'][:]
				v = data.variables['v'][:]
				ukt = units.Quantity(u.values, 'm/s').to(units.knots)
				vkt = units.Quantity(v.values, 'm/s').to(units.knots)

				ws = mpcalc.wind_speed(data['u'], data['v']).metpy.convert_units('knots')

				t = data.variables['t'][:]
				t_c = units.Quantity(t.values, 'K').to(units.celsius)

				z = data.variables['gh'][:]
				z = units.Quantity(z.values, 'm^2/s^2')
				gpm = mpcalc.geopotential_to_height(z)

				gpm_center = round((gpm.magnitude.min() + gpm.magnitude.max())/2, 10)
				gpm_levels = np.arange(gpm_center - 500, gpm_center + 500, 36)

				t_center = round((t_c.magnitude.min() + t_c.magnitude.max())/2, 10)
				t_levels = np.arange(t_center - 20, t_center + 20, 2)

				title_desc = f'Wind Speed (kts, shaded),\nTemperature (°C, colored contours), and Geopotential Height (m, contours)'
				title_info = f'{level} hPa\nValid at {np.datetime_as_string(data["valid_time"].values, timezone="UTC")[:-11].replace("T", " ")}'

			crs = ccrs.PlateCarree(central_longitude = -60)

			smoothing = False
			zoom_fact = 0 #2=zoomed

			center = (53.895496, -166.537517)
			rad = 0.75

			extent = [center[1]+rad, center[1]-rad, center[0]+rad, center[0]-rad]

			#extent = [-104.3+zoom_fact, -100.3-zoom_fact, 28+zoom_fact, 31-zoom_fact] #WESN

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
			gl.xlabels_bottom = True
			gl.ylabels_left = True
			gl.ylabels_right= False
			gl.xlines = True
			gl.ylines = True
			gl.xlocator = mticker.FixedLocator(np.arange(extent[0], extent[1], 1))
			gl.ylocator = mticker.FixedLocator(np.arange(extent[2], extent[3], 1))
				
			# levels = np.arange(0,round(math.ceil(ws.magnitude.max()/10))*10,2)
			levels = np.arange(0, 70, 2)

			#levels = np.arange(0, 150, 2)

			# if float(sel_level) > 500:
			# 	levels = np.arange(-20, 100, 3)
			# else:
			# 	levels = np.arange(-120, 20, 3)



			if smoothing == True:
				wind_contours = ax.contourf(lon, lat, 
				                            mpcalc.smooth_n_point(ws, 9, 1), 
				                            levels = levels, 
				                            transform = ccrs.PlateCarree(), 
				                            cmap=cm.jet, 
				                            extend='max',
				                            zorder=0)
			else:
				wind_contours = ax.contourf(lon, lat, 
				                            ws, 
				                            levels = levels, 
				                            transform = ccrs.PlateCarree(), 
				                            cmap=cm.jet, 
				                            extend='max',
				                            zorder=0)
			cb = plt.colorbar(wind_contours, 
			                  ax=ax, 
			                  ticks = levels[::2], 
			                  orientation="vertical", 
			                  shrink=0.77)
			                  #extendrect=True)
			cb.set_label('Wind Gust Speed (kts)', size='x-large')

			# contours = plt.contour(lon, lat, gpm, levels=gpm_levels, colors="gray", alpha=1, 
			# 	transform=ccrs.PlateCarree(), zorder=11, linewidths=1.8)
			# plt.clabel(contours, inline=1, fontsize=10, fmt="%i", zorder=6)

			# t_contours = plt.contour(lon, lat, t_c, levels=t_levels, cmap=cm.RdYlBu_r, alpha=0.85, 
			# 	transform=ccrs.PlateCarree(), zorder=13, linestyles='dashed', linewidths=1.2,)
			# plt.clabel(t_contours, inline=1, fontsize=10, fmt="%i", zorder=6)

			if wind_display == "Vectors":
				if idx > 0:
					obj = ax.quiver(lon[::idx, ::idx], lat[::idx, ::idx], ukt[::idx, ::idx], vkt[::idx, ::idx], color='greenyellow', scale=600, transform=ccrs.PlateCarree())
				else:
					obj = ax.quiver(lon, lat, ukt, vkt, color='greenyellow', transform=ccrs.PlateCarree())

				obj.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='black')])
			elif wind_display == "Barbs":
				if idx > 0:
					obj = ax.barbs(lon[::idx, ::idx], lat[::idx, ::idx], ukt[::idx, ::idx], vkt[::idx, ::idx], color='greenyellow', length=6, transform=ccrs.PlateCarree())
				else:
					obj = ax.barbs(lon, lat, ukt, vkt, color='greenyellow', length=6, transform=ccrs.PlateCarree())
				obj.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='black')])
			else:
				pass

			airport = (-166.537517, 53.895496)
			plt.plot(airport[0], airport[1], marker='*', markersize=20, markerfacecolor='greenyellow', markeredgecolor='black', transform=ccrs.PlateCarree(), zorder=15)
			# incident_label = plt.annotate("Incident Site", (airport[0], airport[1]-0.2), color='greenyellow', fontweight='bold', horizontalalignment='center', transform=ccrs.PlateCarree(), zorder=15)
			# incident_label.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='black')])

			plot_towns(ax,
			           (extent[0], extent[1]),
			           (extent[2], extent[3]),
			           scale_rank=12)

			# titlestr = f"{title_desc}\n{title_info}"
			# plt.title(titlestr)#, {'fontsize': 32})
			plt.title(f'HRRR Reanalysis {title_desc}', loc='left', fontweight='bold', fontsize=15)
			plt.title(title_info, loc='right')#, {'fontsize': 32})

			draw_logo(ax)

			day_dir = sel_time[:10].replace("-", "_")

			time_dir = sel_time[11:16].replace(":", "") + "UTC"

			#output_dir = os.path.join(save_dir, day_dir, time_dir)

			output_dir = os.path.join(save_dir, wind_display, sel_level)

			if not os.path.exists(output_dir):
				os.makedirs(output_dir)

			#file_name = "PlanView_" + PRODUCT + "_" + sel_level[:-2] + "mb_" + time_dir +"_ERA5.png"
			file_name = "PlanView_" + day_dir + "_" + time_dir + "_" + wind_display + "_HRRR.png"

			plt.savefig(os.path.join(output_dir, file_name), bbox_inches="tight", dpi=200)

			plt.close()

			print("				Saved!")

print("DONE!!!!!!!!!!!!!!!!! FINALLY!!!!!!!!!!!!!!")



