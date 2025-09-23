import os
import cdsapi


c = cdsapi.Client()

data_path = '/Users/rpurciel/Documents/data'
year = '2019'
month = '03'
days = ['02', '03']
times = ['00:00', '01:00', '02:00',
        '03:00', '04:00', '05:00',
        '06:00', '07:00', '08:00',
        '09:00', '10:00', '11:00',
        '12:00', '13:00', '14:00',
        '15:00', '16:00', '17:00',
        '18:00', '19:00', '20:00',
        '21:00', '22:00', '23:00']

bbox = [12, 25, -12, 50]

for day in days:
    for time in times:

        ua_fname = f'{year}_{month}_{day}_{time.replace(":", "")}UTC_UA_ERA5.grib'
        sfc_fname = f'{year}_{month}_{day}_{time.replace(":", "")}UTC_SFC_ERA5.grib'
        land_fname = f'{year}_{month}_{day}_{time.replace(":", "")}UTC_ERA5Land.grib'

        ua_path = os.path.join(data_path, ua_fname)
        sfc_path = os.path.join(data_path, sfc_fname)
        land_path = os.path.join(data_path, 'ERA5Land', land_fname)

        # c.retrieve(
        #     'reanalysis-era5-pressure-levels',
        #     {
        #         'product_type': 'reanalysis',
        #         'format': 'grib',
        #         'variable': [
        #             'divergence', 'fraction_of_cloud_cover', 'geopotential',
        #             'ozone_mass_mixing_ratio', 'potential_vorticity', 'relative_humidity',
        #             'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content', 'specific_humidity',
        #             'specific_rain_water_content', 'specific_snow_water_content', 'temperature',
        #             'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity',
        #             'vorticity',
        #         ],
        #         'pressure_level': [
        #             '1', '2', '3',
        #             '5', '7', '10',
        #             '20', '30', '50',
        #             '70', '100', '125',
        #             '150', '175', '200',
        #             '225', '250', '300',
        #             '350', '400', '450',
        #             '500', '550', '600',
        #             '650', '700', '750',
        #             '775', '800', '825',
        #             '850', '875', '900',
        #             '925', '950', '975',
        #             '1000',
        #         ],
        #         'year': year,
        #         'month': month,
        #         'day': day,
        #         'time': time,
        #         'area': bbox
        #     },
        #     ua_path)

        # c.retrieve(
        #     'reanalysis-era5-single-levels',
        #     {
        #         'product_type': 'reanalysis',
        #         'format': 'grib',
        #         'variable': [
        #             '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
        #             '2m_temperature', 'geopotential', 'surface_pressure', '10m_wind_gust_since_previous_post_processing',
        #             'instantaneous_10m_wind_gust'
        #         ],
        #         'year': year,
        #         'month': month,
        #         'day': day,
        #         'time': time,
        #         'area': bbox
        #     },
        #     sfc_path)

        c.retrieve(
            'reanalysis-era5-land',
            {
                'format': 'grib',
                'variable': [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
                    '2m_temperature', 'geopotential', 'surface_pressure',
                ],
                'year': year,
                'month': month,
                'day': day,
                'time': time,
                'area': bbox
            },
            land_path)



