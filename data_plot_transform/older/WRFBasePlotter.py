
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
from netCDF4 import Dataset
import wrf
import pint

#STD_LEVELS = [1000, 950, 900, 850, 700, 500, 300, 250]
STD_LEVELS = [700, 725, 750, 775, 800]
#STD_LEVELS = ['SFC']

STD_HOURS = []

PRODUCT = "VerticalVelocity"

idx = 10

wind_display = "None"

warnings.filterwarnings("ignore")

save_dir = "/Users/rpurciel/Documents/case/WRF/VV Plan View/"

path = "/Users/rpurciel/Documents/case/WRF/Data/"

files = sorted(glob.glob(path + "wrfout_d04*"))

print("Starting...")

for file in files:

    sfc_file = False

    print("Selected File ", file)

    raw_data = []
    raw_data += [xr.open_dataset(file, engine="netcdf4")]
    raw_data_nc4 = Dataset(file)
    raw_data = xr.merge(raw_data)

    num_times = 1

    if sfc_file:
        num_levels = 1
        sel_levels = ["SFC"]
    else:
        num_levels = len(raw_data.bottom_top)
        sel_levels = STD_LEVELS

    for time in range(0, num_times, 1):

        timesel_data = raw_data

        sel_time = str(timesel_data['XTIME'].values[0])

        print("        Time: ", sel_time)

        for level in sel_levels:

            if sfc_file:
                levelsel_data = timesel_data
                sel_level = 'SFCXX'

            else:
                levelsel_data = timesel_data
                sel_level = str(level)

            print("            Level: ", sel_level)

            data = levelsel_data

            #print(data.variables)

            lat = data.XLAT.values[0]
            lon = data.XLONG.values[0]

            if sfc_file:

                u = units.Quantity(data.U10.values, 'm/s')
                v = units.Quantity(data.V10.values, 'm/s')
                ws_ms = mpcalc.wind_speed(u, v)
                ws = ws_ms.to(units.knots)

                ws = ws[0]

                temp_k = units.Quantity(data.T2.values, 'K')

                ukt = u.to(units.knots)
                vkt = v.to(units.knots)

                z = units.Quantity(data.PHP.values, 'm^2/s^2')

                gpm = mpcalc.geopotential_to_height(z[0])

                t_c = temp_k.to(units.celsius)
                t_c = t_c[0]

            else:

                ureg = pint.UnitRegistry()

                # # p_pa = units.Quantity(data.P.data + data.PB.data, 'pascal')
                # # p = p_pa.to(ureg('hectopascal'))

                # # u_stg = data.U

                # # u_dstg = wrf.destagger(data.U, 2, meta=True)
                # # v_dtsg = wrf.destagger(data.V, 2, meta=True)

                # # print(p, u_stg, u_dstg)
                # # print(p.shape, u_stg.shape, u_dstg.shape)

                # lat = 41.361298
                # lon = -111.759610

                # latitude = wrf.getvar(raw_data_nc4, "lat")
                # longitude = wrf.getvar(raw_data_nc4, "lon")
                
                # # Find the nearest neightbor from the maximum absolute differences...
                # abslat = np.abs(latitude-lat)
                # abslon = np.abs(longitude-lon)

                # c = np.maximum(abslon, abslat)
                
                # # Get the index of the nearest point
                # # (Again, the values are backwards you might expect because the 
                # # coordinates are in (y, x) order.)
                # ([idx_y], [idx_x]) = np.where(c == np.min(c))
                
                # #Now that we have the index location of the nearest point to the requested lat/lon, 
                # #we can use the select function sel to get all the data in the height coordinate at a single point in the dataset.
                
                # point_data = data.sel(south_north=idx_y, west_east=idx_x)

                # col_press_pa = units.Quantity(point_data.P.data + point_data.PB.data, 'pascal')
                # col_press_hpa = col_press_pa.to(ureg('hectopascal'))

                # print(col_press_hpa.values)


                vv_da = wrf.interplevel(wrf.getvar(raw_data_nc4, "wa", units='ft s-1'), wrf.getvar(raw_data_nc4, "pressure"), level)
                vv = units.Quantity(vv_da.values, "ft/s")

                ws_da, _ = wrf.interplevel(wrf.getvar(raw_data_nc4, "uvmet_wspd_wdir", units='kt'), wrf.getvar(raw_data_nc4, "pressure"), level)
                ws = units.Quantity(ws_da.values, 'knots')

                ukt_da, vkt_da = wrf.interplevel(wrf.getvar(raw_data_nc4, "uvmet", units='kt'), wrf.getvar(raw_data_nc4, "pressure"), level)
                ukt = units.Quantity(ukt_da.values, 'knots')
                vkt = units.Quantity(vkt_da.values, 'knots')

                t_c_da = wrf.interplevel(wrf.getvar(raw_data_nc4, "tc"), wrf.getvar(raw_data_nc4, "pressure"), level)
                t_c = units.Quantity(t_c_da.values, 'celsius')

                gpm_da = wrf.interplevel(wrf.getvar(raw_data_nc4, "z", units='m'), wrf.getvar(raw_data_nc4, "pressure"), level)
                gpm = units.Quantity(gpm_da.values, 'meters')

            crs = ccrs.PlateCarree(central_longitude = -60)

            smoothing = False
            #extent = [110, 155, -10, -45] #WESN, australia

            # zoom_fact = 4 #2=zoomed

            # extent = [30.5+zoom_fact, 42-zoom_fact, -1.5+zoom_fact, 8.5-zoom_fact]

            extent = [-112.378, -111.313, 40.954, 41.662] #Cessna Eden UT

            #Bell Kenya

            # extent = [37.149, 34.933, 2.395, 4.608] #d03

            # extent = [36.374, 35.708, 3.169, 3.835] #d04

            # extent = [36.140, 35.942, 3.403, 3.601] #d05

            fig = plt.figure(figsize=(22,16))
            ax = plt.axes(projection = ccrs.PlateCarree())
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            states = NaturalEarthFeature(category="cultural", scale="50m",
                                                  facecolor="none",
                                                  name="admin_1_states_provinces")
            ax.add_feature(states, linewidth=1.0, edgecolor="black")
            ax.coastlines('50m', linewidth=1.5)
            ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), linestyle='-', linewidth=0.5, alpha=1,edgecolor='white',facecolor='none')
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
                
            #levels = np.arange(0,(round(math.ceil(ws.max()/10))*10)+2,2)

            levels = np.arange(0, 60, 2)

            #levels = np.arange(0, 150, 2)

            # print(np.nan_to_num(vv.magnitude, copy=False).min(), np.nan_to_num(vv.magnitude, copy=False).max())
            # vv_levels = np.arange(np.nan_to_num(vv.magnitude, copy=False).min(), np.nan_to_num(vv.magnitude, copy=False).max(), 0.5)
            vv_levels = np.arange(-40, 40, 4)

            # if float(sel_level) > 500:
            
            #     levels = np.arange(-20, 100, 3)
            # else:
            #     levels = np.arange(-120, 20, 3)

            # gpm_center = round((gpm.magnitude.min() + gpm.magnitude.max()/2), 10)

            # print(gpm.magnitude)

            # gpm_levels = np.arange(gpm_center - 500, gpm_center + 500, 36)

            # t_center = round((t_c.magnitude.min() + t_c.magnitude.max())/2, 10)

            # t_levels = np.arange(t_center - 20, t_center + 20, 2)


            if smoothing == True:
                # wind_contours = ax.contourf(lon, lat, mpcalc.smooth_n_point(ws, 9, 1), levels = levels, transform = ccrs.PlateCarree(), cmap=cm.jet, zorder=0)
                pass
            else:
                # wind_contours = ax.contourf(lon, lat, ws.magnitude, levels = levels, transform = ccrs.PlateCarree(), cmap=cm.jet, zorder=0)
                wind_contours = ax.contourf(lon, lat, vv.magnitude, levels=vv_levels, transform=ccrs.PlateCarree(), cmap=cm.bwr, zorder=0)
            cb = plt.colorbar(wind_contours, ax=ax, ticks = vv_levels[::2], orientation="vertical", shrink=0.77, extend="both", extendrect=True)
            # cb.set_label('Wind Speed (kts)', size='x-large')
            cb.set_label('Vertical Velocty (ft/sec)', size='x-large')

            # contours = plt.contour(lon, lat, gpm, levels=gpm_levels, colors="gray", alpha=1, 
            #     transform=ccrs.PlateCarree(), zorder=11, linewidths=1.8)
            # plt.clabel(contours, inline=1, fontsize=10, fmt="%i", zorder=6)

            # t_contours = plt.contour(lon, lat, t_c.magnitude, levels=t_levels, cmap=cm.RdYlBu_r, alpha=0.85, 
            #     transform=ccrs.PlateCarree(), zorder=13, linestyles='dashed', linewidths=1.2,)
            # plt.clabel(t_contours, inline=1, fontsize=10, fmt="%i", zorder=6)

            # ukt = ukt[0, :, :]
            # vkt = vkt[0, :, :]

            print(lon.ndim, lat.ndim, ukt.magnitude.ndim, vkt.magnitude.ndim)
            print(lon.shape, lat.shape, ukt.magnitude.shape, vkt.magnitude[::8].shape)
            if wind_display == "Vectors":
                if idx > 0:
                    obj = ax.quiver(lon[::idx], lat[::idx], ukt.magnitude[::idx], vkt.magnitude[::idx], color='greenyellow', transform=ccrs.PlateCarree())
                else:
                    obj = ax.quiver(lon, lat, ukt, vkt, color='greenyellow', transform=ccrs.PlateCarree())

                obj.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='black')])
            elif wind_display == "Barbs":
                if idx > 0:
                    obj = ax.barbs(lon[::idx], lat[::idx], ukt.magnitude[::idx], vkt.magnitude[::idx], color='greenyellow', length=6, transform=ccrs.PlateCarree())
                else:
                    obj = ax.barbs(lon, lat, ukt, vkt, color='greenyellow', length=6, transform=ccrs.PlateCarree())
                obj.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='black')])
            else:
                pass

            if sfc_file:
                titlestr = "Wind Speed (kts, shaded) and Temperature (°C, colored contours)\n2 meters above Surface, Valid at " + np.datetime_as_string(data['XTIME'].values[0], timezone="UTC")[:-11].replace("T", " ")
            else:
                titlestr = "Wind Speed (kts, shaded) and Temperature (°C, colored contours)\n isobaricInhPa (" + str(level) +" hPa), Valid at " + np.datetime_as_string(data['XTIME'].values[0], timezone="UTC")[:-11].replace("T", " ")

            titlestr = "Vertical Velocity (ft/sec, shaded)\n isobaricInhPa (" + str(level) +" hPa), Valid at " + np.datetime_as_string(data['XTIME'].values[0], timezone="UTC")[:-11].replace("T", " ")            #titlestr = "Wind Speed (kts, shaded), Temperature (°C, colored contours), and Geopotential Height (m, contours)\n" + t.attrs['GRIB_typeOfLevel'] + " (" + str(data['isobaricInhPa'].values) +" hPa), Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")[:-11].replace("T", " ")
            #titlestr = "Wind Speed (kts, shaded) and Temperature (°C, colored contours)\n" + t.attrs['GRIB_typeOfLevel'] + " (" + str(data['isobaricInhPa'].values) +" hPa), Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")[:-11].replace("T", " ")
            # titlestr = "Wind Speed (kts, shaded) and Temperature (°C, colored contours)\n2 meters above Surface, Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")[:-11].replace("T", " ")
            #titlestr = "Temperature (°C) and 500 mb Geopotential Height (m)\n" + t.attrs['GRIB_typeOfLevel'] + " (" + str(data['isobaricInhPa'].values) +" hPa), Valid at " + np.datetime_as_string(data['valid_time'].values, timezone="UTC")
            #titlestr = "Winds and Wind Speed\n" + data.attrs['RDA_DATASET_GROUP'] + ", Valid at " + np.datetime_as_string(data['time'].values, timezone="UTC")
            plt.title(titlestr)#, {'fontsize': 32})

            incident = (-111.759610, 41.361298)
            plt.plot(incident[0], incident[1], marker='X', markersize=14, markerfacecolor='greenyellow', markeredgecolor='black', transform=ccrs.PlateCarree(), zorder=15)
            incident_label = plt.annotate("Accident Site", (incident[0], incident[1]-0.03), color='greenyellow', fontweight='bold', horizontalalignment='center', transform=ccrs.PlateCarree(), zorder=15)
            incident_label.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='black')])

            # incident = (36.040758, 3.501691)
            # plt.plot(incident[0], incident[1], marker='X', markersize=14, markerfacecolor='greenyellow', markeredgecolor='black', transform=ccrs.PlateCarree(), zorder=15)
            # incident_label = plt.annotate("Incident Site", (incident[0], incident[1]-0.007), color='greenyellow', fontweight='bold', horizontalalignment='center', transform=ccrs.PlateCarree(), zorder=15)
            # incident_label.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='black')])

            # airport = (35.605500, 3.121481)
            # plt.plot(airport[0], airport[1], marker='X', markersize=14, markerfacecolor='greenyellow', markeredgecolor='black', transform=ccrs.PlateCarree(), zorder=15)
            # airport_label = plt.annotate("HKLO", (airport[0], airport[1]-0.007), color='greenyellow', fontweight='bold', horizontalalignment='center', transform=ccrs.PlateCarree(), zorder=15)
            # airport_label.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='black')])

            day_dir = sel_time[:10].replace("-", "_")

            time_dir = sel_time[11:16].replace(":", "") + "UTC"

            #output_dir = os.path.join(save_dir, day_dir, time_dir)

            output_dir = os.path.join(save_dir, wind_display, sel_level)

            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            #file_name = "PlanView_" + PRODUCT + "_" + sel_level[:-2] + "mb_" + time_dir +"_ERA5.png"
            file_name = "PlanView_" + PRODUCT + "_" + day_dir + "_" + time_dir +"_WRF.png"

            plt.savefig(os.path.join(output_dir, file_name), bbox_inches="tight", dpi=200)

            plt.close()

            print("                Saved!")

print("DONE!!!!!!!!!!!!!!!!! FINALLY!!!!!!!!!!!!!!")



