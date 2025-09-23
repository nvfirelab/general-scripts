
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

from __internal_funcs import plot_towns, make_highly_visible, draw_logo

#prod = ['temp_c', 'wind_cf', 'wind_display']
prod = ['wind_cf', 'wind_display']

STD_LEVELS = [1000, 950, 900, 850, 700, 500, 300, 250]
#STD_LEVELS = ['SFC']

STD_HOURS = []

PRODUCT = "Wind"

idx = 2

wind_display = "Vectors"

warnings.filterwarnings("ignore")

save_dir = "/Users/rpurciel/Documents/case/ERA5 Winds/"

path = "/Users/rpurciel/Documents/case/ERA5/data/"

files = sorted(glob.glob(path + "*UA_ERA5.grib"))

print("Starting...")

for file in files:

	if ("SFC" in file) or ("Land" in file):
		sfc_file = True
	else:
		sfc_file = False

	print("Selected File ", file)

	raw_data = []
	raw_data += [xr.open_dataset(file, engine="cfgrib")]
	#raw_data += [xr.open_dataset(file, engine="cfgrib", filter_by_keys={'short_name': 'i10fg'})]
	raw_data = xr.merge(raw_data, compat='override')

	num_times = 1

	if sfc_file:
		num_levels = 1
		sel_levels = ["SFC"]
	else:
		num_levels = len(raw_data.isobaricInhPa)
		sel_levels = STD_LEVELS

	for time in range(0, num_times, 1):

		timesel_data = raw_data

		sel_time = str(timesel_data['time'].values)

		print("		Time: ", sel_time)

		for level in sel_levels:

			if sfc_file:
				levelsel_data = timesel_data
				sel_level = 'SFCXX'

			else:
				levelsel_data = timesel_data.sel(isobaricInhPa=level)
				sel_level = levelsel_data['isobaricInhPa'].values

			print("			Level: ", sel_level)

			# extent = [250,359,0,65]

			data = levelsel_data

			# nd1 = levelsel_data.where(levelsel_data.latitude <= 65, drop=True) 
			# nd2 = nd1.where(levelsel_data.latitude >= -0, drop=True)
			# nd3 = nd2.where(levelsel_data.longitude <= 359, drop=True)
			# nd4 = nd3.where(levelsel_data.longitude >= 250, drop=True)

			# data = nd4

			#print(data.variables)

			lat = data.variables['latitude'][:]
			lon = data.variables['longitude'][:]

			if sfc_file:
				u = data.variables['u10'][:]
				v = data.variables['v10'][:]
				ws = mpcalc.wind_speed(data['u10'], data['v10']).metpy.convert_units('knots')
				#ws = units.Quantity(data['i10fg'].values, 'm/s').to(units.knots)
				t = data.variables['t2m'][:]
				#rh = data.variables
			else:
				u = data.variables['u'][:]
				v = data.variables['v'][:]
				ws = mpcalc.wind_speed(data['u'], data['v']).metpy.convert_units('knots')
				t = data.variables['t'][:]
				rh = data.variables['r'][:]
				cc = data.variables['cc'][:]

			ukt = units.Quantity(u.values, 'm/s').to(units.knots)
			vkt = units.Quantity(v.values, 'm/s').to(units.knots)

			# z = data.variables['z'][:]

			# z = units.Quantity(z.values, 'm^2/s^2')

			# gpm = mpcalc.geopotential_to_height(z)

			t_c = units.Quantity(t.values, 'K').to(units.celsius)

			crs = ccrs.PlateCarree(central_longitude = -60)

			smoothing = False
			#extent = [110, 155, -10, -45] #WESN, australia

			zoom_fact = 2 #2=zoomed

			extent = [30.5+zoom_fact, 42-zoom_fact, -1.5+zoom_fact, 8.5-zoom_fact] #BK
			#extent = [34.75+zoom_fact, 39.25-zoom_fact, -5.5+zoom_fact, -1.25-zoom_fact] #Mt Kilimanjaro.

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
			gl.ylocator = mticker.FixedLocator(np.arange(extent[2], extent[3], 1))
			    
			#levels = np.arange(0,(round(math.ceil(ws.max()/10))*10)+2,2)

			#levels = np.arange(0, 150, 2)

			# if float(sel_level) > 500:
			
			# 	levels = np.arange(-20, 100, 3)
			# else:
			# 	levels = np.arange(-120, 20, 3)

			# gpm_center = round((gpm.magnitude.min() + gpm.magnitude.max()/2), 10)

			# gpm_levels = np.arange(gpm_center - 500, gpm_center + 500, 36)

			# t_center = round((np.nanmin(t_c.magnitude) + np.nanmax(t_c.magnitude))/2, 10)

			# t_levels = np.arange(t_center - 20, t_center + 20, 2)

			prodstr = ""
			prodabvr = ""

			if "cldcvr_cf" in prod:
				levels = np.arange(0, 1.1, 0.1)
				if smoothing == True:
				    rh_contours = ax.contourf(lon, lat, 
				                                mpcalc.smooth_n_point(cc, 9, 1), 
				                                levels = levels, 
				                                transform = ccrs.PlateCarree(), 
				                                cmap=cm.Blues, 
				                                zorder=0)
				else:
				    rh_contours = ax.contourf(lon, lat, 
				                                cc, 
				                                levels = levels, 
				                                transform = ccrs.PlateCarree(), 
				                                cmap=cm.Blues, 
				                                vmin=levels.min(),
				                                vmax=levels.max(),
				                                zorder=0)
				cb = plt.colorbar(rh_contours, 
				                  ax=ax, 
				                  ticks = levels[::2], 
				                  orientation="vertical", 
				                  shrink=0.77,)
				cb.set_label('Fraction of Cloud Cover', size='x-large')

				prodstr += "#Fraction of Cloud Cover (shaded)"
				prodabvr += "_CC"

			if "rh_cf" in prod:
				levels = np.arange(0, 105, 5)
				if smoothing == True:
				    rh_contours = ax.contourf(lon, lat, 
				                                mpcalc.smooth_n_point(rh, 9, 1), 
				                                levels = levels, 
				                                transform = ccrs.PlateCarree(), 
				                                cmap=cm.terrain_r, 
				                                extend="max",
				                                zorder=0)
				else:
				    rh_contours = ax.contourf(lon, lat, 
				                                rh, 
				                                levels = levels, 
				                                transform = ccrs.PlateCarree(), 
				                                cmap=cm.terrain_r, 
				                                vmin=levels.min(),
				                                vmax=levels.max(),
				                                extend="max",
				                                zorder=0)
				cb = plt.colorbar(rh_contours, 
				                  ax=ax, 
				                  ticks = levels[::2], 
				                  orientation="vertical", 
				                  shrink=0.77,
				                  extendrect=True,)
				cb.set_label('Relative Humidty (%)', size='x-large')

				prodstr += "#Relative Humidity (%, shaded)"
				prodabvr += "_RH"

			if "wind_cf" in prod:
				levels = np.arange(0, 36, 2)
				if smoothing == True:
				    wind_contours = ax.contourf(lon, lat, 
				                                mpcalc.smooth_n_point(ws, 9, 1), 
				                                levels = levels, 
				                                transform = ccrs.PlateCarree(), 
				                                extend='max', 
				                                cmap=cm.jet, 
				                                zorder=0)
				else:
				    wind_contours = ax.contourf(lon, lat, 
				                                ws, 
				                                levels = levels, 
				                                transform = ccrs.PlateCarree(), 
				                                cmap=cm.jet, 
				                                vmin=levels.min(),
				                                vmax=levels.max(),
				                                extend='max', 
				                                zorder=0)
				cb = plt.colorbar(wind_contours, 
				                  ax=ax, 
				                  ticks = levels[::2], 
				                  orientation="vertical", 
				                  shrink=0.77, 
				                  extend="max")
				cb.set_label('Wind Speed (kts)', size='x-large')

				prodstr += "#Wind Speed (kts, shaded)"
				prodabvr += "_WS"

			if "vv_cf" in prod:
				levels = np.arange(0, 62, 2)
				if smoothing == True:
				    wind_contours = ax.contourf(lon, lat, 
				                                mpcalc.smooth_n_point(ws, 9, 1), 
				                                levels = levels, 
				                                transform = ccrs.PlateCarree(), 
				                                extend='max', 
				                                cmap=cm.jet, 
				                                zorder=0)
				else:
				    wind_contours = ax.contourf(lon, lat, 
				                                ws, 
				                                levels = levels, 
				                                transform = ccrs.PlateCarree(), 
				                                cmap=cm.jet, 
				                                vmin=levels.min(),
				                                vmax=levels.max(),
				                                extend='max', 
				                                zorder=0)
				cb = plt.colorbar(wind_contours, 
				                  ax=ax, 
				                  ticks = levels[::2], 
				                  orientation="vertical", 
				                  shrink=0.77, 
				                  extend="max")
				cb.set_label('Wind Speed (kts)', size='x-large')

				prodstr += "#Wind Speed (kts, shaded)"
				prodabvr += "_WS"

			if "gpm_c" in prod:

				gpm_center = round((gpm.magnitude.min() + gpm.magnitude.max()/2), 10)
				gpm_levels = np.arange(gpm_center - 500, gpm_center + 500, 36)

				contours = plt.contour(lon, lat, 
				                       gpm_levels, 
				                       levels=levels, 
				                       colors="gray", 
				                       alpha=1, 
									   transform=ccrs.PlateCarree(), 
									   zorder=11, 
									   linewidths=1.8)
				plt.clabel(contours, inline=1, fontsize=10, fmt="%i", zorder=6)

				prodstr += "#Geopotential Height (m, contours)"
				prodabvr += "_GPM"


			if "temp_c" in prod:

				t_min = -30 if sel_level > 500 else -70
				t_max = 45 if sel_level > 500 else 20
				t_levels = np.arange(t_min, t_max, 1)

				t_contours = plt.contour(lon, lat, 
				                         t_c, 
				                         levels=t_levels, 
				                         cmap=cm.RdYlBu_r, 
				                         alpha=0.85, 
										 transform=ccrs.PlateCarree(), 
										 zorder=13, 
										 linestyles='dashed', 
										 linewidths=1.2,)
				plt.clabel(t_contours, inline=1, fontsize=10, fmt="%i", zorder=6)

				prodstr += "#Temperature (°C, colored contours)"
				prodabvr += "_T"


			if "wind_display" in prod:
				if wind_display == "Vectors":
					if idx > 0:
						obj = ax.quiver(lon[::idx], lat[::idx], ukt[::idx, ::idx], vkt[::idx, ::idx], color='greenyellow', transform=ccrs.PlateCarree())
					else:
						obj = ax.quiver(lon, lat, ukt, vkt, color='greenyellow', transform=ccrs.PlateCarree())

					obj.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='black')])
					prodabvr += "_Wv"
				elif wind_display == "Barbs":
					if idx > 0:
						obj = ax.barbs(lon[::idx], lat[::idx], ukt[::idx, ::idx], vkt[::idx, ::idx], color='greenyellow', length=6, transform=ccrs.PlateCarree())
					else:
						obj = ax.barbs(lon, lat, ukt, vkt, color='greenyellow', length=6, transform=ccrs.PlateCarree())
					obj.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='black')])
					prodabvr += "_Wb"
				elif wind_display == "Values":
					lat_val, lon_val = np.meshgrid(lat.values, lon.values)
					if idx > 0:
						for ws_row, lat_row, lon_row in zip(ws, lat_val, lon_val):
							for this_ws, this_lat, this_lon in zip(ws_row[::idx], lat_row[::idx], lon_row[::idx]):
								obj = ax.annotate(str(round(this_ws.magnitude, 1)), (this_lon, this_lat), horizontalalignment='right', verticalalignment='top', color='white', clip_box=ax.bbox, fontsize=8, transform=ccrs.PlateCarree(), annotation_clip=False, zorder=30)
								# obj.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='black')])
					else:
						for ws_row, lat_row, lon_row in zip(ws, lat_val, lon_val):
							for this_ws, this_lat, this_lon in zip(ws_row, lat_row, lon_row,):
								obj = ax.annotate(str(round(this_ws.magnitude, 1)), (this_lon, this_lat), horizontalalignment='right', verticalalignment='top', color='white', clip_box=ax.bbox, fontsize=8, transform=ccrs.PlateCarree(), annotation_clip=False, zorder=30)
								obj.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='black')])

					prodabvr += "_Wd"
				else:
					pass

			if sfc_file:
				titlestr = "10m Instantaneous Wind Gust Speed (kts, shaded) and 2m Temperature (°C, colored contours)"
				descstr = "Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")[:-11].replace("T", " ")
			else:
				prod_titles = prodstr.split('#')
				if len(prod_titles[1:]) > 1:
					prod_titles[-1] = "and " + prod_titles[-1]
					if len(prodstr) >= 40:
						prod_titles[round(len(prod_titles)/2)] = '\n' + prod_titles[round(len(prod_titles)/2)]
				titlestr = ", ".join(prod_titles)
				titlestr = titlestr[2:]
				descstr = t.attrs['GRIB_typeOfLevel'] + " (" + str(data['isobaricInhPa'].values) +" hPa), Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")[:-11].replace("T", " ")

			plot_towns(ax,
			           (extent[0], extent[1]),
			           (extent[2], extent[3]),
			           scale_rank=14)

			#titlestr = "Wind Speed (kts, shaded), Temperature (°C, colored contours), and Geopotential Height (m, contours)\n" + t.attrs['GRIB_typeOfLevel'] + " (" + str(data['isobaricInhPa'].values) +" hPa), Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")[:-11].replace("T", " ")
			#titlestr = "Wind Speed (kts, shaded) and Temperature (°C, colored contours)\n" + t.attrs['GRIB_typeOfLevel'] + " (" + str(data['isobaricInhPa'].values) +" hPa), Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")[:-11].replace("T", " ")
			# titlestr = "Wind Speed (kts, shaded) and Temperature (°C, colored contours)\n2 meters above Surface, Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")[:-11].replace("T", " ")
			#titlestr = "Temperature (°C) and 500 mb Geopotential Height (m)\n" + t.attrs['GRIB_typeOfLevel'] + " (" + str(data['isobaricInhPa'].values) +" hPa), Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")
			#titlestr = "Winds and Wind Speed\n" + data.attrs['RDA_DATASET_GROUP'] + ", Valid at " + np.datetime_as_string(data['time'].values, timezone="UTC")
			plt.title(f'ERA5 Reanalysis {titlestr}', loc='left', fontweight='bold', fontsize=15)
			plt.title(descstr, loc='right')#, {'fontsize': 32})

			airport = (36.040758, 3.501691) #BK
			plt.plot(airport[0], airport[1], marker='*', markersize=20, markerfacecolor='greenyellow', markeredgecolor='black', transform=ccrs.PlateCarree(), zorder=15)
			# incident_label = plt.annotate("Incident Site", (airport[0], airport[1]-0.2), color='greenyellow', fontweight='bold', horizontalalignment='center', transform=ccrs.PlateCarree(), zorder=15)
			# incident_label.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='black')])

			draw_logo(ax)

			day_dir = sel_time[:10].replace("-", "_")

			time_dir = sel_time[11:16].replace(":", "") + "UTC"

			#output_dir = os.path.join(save_dir, day_dir, time_dir)

			output_dir = os.path.join(save_dir, str(int(sel_level)))

			if not os.path.exists(output_dir):
				os.makedirs(output_dir)

			#file_name = "PlanView_" + PRODUCT + "_" + sel_level[:-2] + "mb_" + time_dir +"_ERA5.png"
			file_name = "PlanView" + prodabvr + "_" + day_dir + "_" + time_dir +"_ERA5.png"

			plt.savefig(os.path.join(output_dir, file_name), bbox_inches="tight", dpi=300)

			plt.close()

			print("				Saved!")

print("DONE!!!!!!!!!!!!!!!!! FINALLY!!!!!!!!!!!!!!")



