#!/usr/bin/env python

"""
Turn the MAESPA input file into a CABLE netcdf file. We can swap MAESPA data
for the raw data later.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (04.08.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import netCDF4 as nc
import datetime

def main(in_fname, out_fname, co2x):

    DEG_2_KELVIN = 273.15
    SW_2_PAR = 2.3
    PAR_2_SW = 1.0 / SW_2_PAR

    df = pd.read_csv(in_fname)

    now = datetime.datetime.now()
    n_time_steps = len(df)
    ndim = 1

    f = nc.Dataset(out_fname, 'w', format='NETCDF4')
    f.description = 'EucFACE met data, created by Martin De Kauwe'
    f.history = "Created %s" % (now)
    f.contact = "mdekauwe@gmail.com"

    # dimensions
    f.createDimension('time', n_time_steps)
    f.createDimension('x', ndim)
    f.createDimension('y', ndim)

    # variables
    time = f.createVariable('time', 'f8', ('time',))
    time.units = "seconds since %s 00:00:00" % (df.Date[0])

    x = f.createVariable('x', 'f8', ('x',))
    x.units = ""

    y = f.createVariable('y', 'f8', ('y',))
    y.units = ""

    latitude = f.createVariable('latitude', 'f8', ('y', 'x',))
    latitude.units = "degrees_north"
    latitude._FillValue = -9999.
    latitude.long_name = "Latitude"

    longitude = f.createVariable('longitude', 'f8', ('y', 'x',))
    longitude.units = "degrees_east"
    longitude._FillValue = -9999.
    longitude.long_name = "Longitude"

    SWdown = f.createVariable('SWdown', 'f8', ('time', 'y', 'x',))
    SWdown.units = "W/m^2"
    SWdown.missing_value = -9999.
    SWdown.long_name = "Surface incident shortwave radiation"
    SWdown.CF_name = "surface_downwelling_shortwave_flux_in_air"

    Tair = f.createVariable('Tair', 'f8', ('time', 'y', 'x',))
    Tair.units = "K"
    Tair.missing_value = -9999.
    Tair.long_name = "Near surface air temperature"
    Tair.CF_name = "surface_temperature"

    Rainf = f.createVariable('Rainf', 'f8', ('time', 'y', 'x',))
    Rainf.units = "mm/s"
    Rainf.missing_value = -9999.
    Rainf.long_name = "Rainfall rate"
    Rainf.CF_name = "precipitation_flux"

    Qair = f.createVariable('Qair', 'f8', ('time', 'y', 'x',))
    Qair.units = "kg/kg"
    Qair.missing_value = -9999.
    Qair.long_name = "Near surface specific humidity"
    Qair.CF_name = "surface_specific_humidity"

    Wind = f.createVariable('Wind', 'f8', ('time', 'y', 'x',))
    Wind.units = "m/s"
    Wind.missing_value = -9999.
    Wind.long_name = "Scalar windspeed" ;
    Wind.CF_name = "wind_speed"

    PSurf = f.createVariable('PSurf', 'f8', ('time', 'y', 'x',))
    PSurf.units = "Pa"
    PSurf.missing_value = -9999.
    PSurf.long_name = "Surface air pressure"
    PSurf.CF_name = "surface_air_pressure"

    LWdown = f.createVariable('LWdown', 'f8', ('time', 'y', 'x',))
    LWdown.units = "W/m^2"
    LWdown.missing_value = -9999.
    LWdown.long_name = "Surface incident longwave radiation"
    LWdown.CF_name = "surface_downwelling_longwave_flux_in_air"

    CO2 = f.createVariable('CO2', 'f8', ('time', 'y', 'x',))
    CO2.units = "umol/mol"
    CO2.missing_value = -9999.
    CO2.long_name = ""
    CO2.CF_name = ""

    elevation = f.createVariable('elevation', 'f8', ('y', 'x',))
    elevation.units = "m" ;
    elevation.missing_value = -9999.
    elevation.long_name = "Site elevation above sea level" ;

    reference_height = f.createVariable('reference_height', 'f8', ('y', 'x',))
    reference_height.units = "m"
    reference_height.missing_value = -9999.
    reference_height.long_name = "Measurement height on flux tower"

    # data
    x = 150.740278 # Ellsworth 2017, NCC
    y = -33.617778 # Ellsworth 2017, NCC
    #time =
    latitude =
    longitude =
    SWdown = df.PAR * PAR_2_SW
    Tair = df.TAIR.values + DEG_2_KELVIN
    Rainf = df.PPT.values
    #Qair =
    Wind = df.WIND.values
    PSurf = df.PRESS.values
    #LWdown =
    if co2x == "amb":
        CO2 = df["Ca.A"].values
    elif co2x == "ele":
        CO2 = df["Ca.E"].values
    elevation = 23.0 # Ellsworth 2017, NCC
    #reference_height =


    f.close()



if __name__ == "__main__":

    in_fname = "raw_data/eucdata.csv"
    for co2_conc in ["amb", "ele"]:
        out_fname = "EucFACE_met_%s.nc" % (co2_conc)
        main(in_fname, out_fname, co2x=co2_conc)
