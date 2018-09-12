#!/usr/bin/env python

"""
Turn the MAESPA input file into a CABLE netcdf file. Aim to swap MAESPA data
for the raw data later when I have more time...

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

def main(in_fname, out_fname, co2_conc):

    DEG_2_KELVIN = 273.15
    SW_2_PAR = 2.3
    PAR_2_SW = 1.0 / SW_2_PAR

    df = pd.read_csv(in_fname)

    ndim = 1
    n_timesteps = len(df)
    times = []
    secs = 0.0
    for i in range(n_timesteps):
        times.append(secs)
        secs += 1800.

    # create file and write global attributes
    f = nc.Dataset(out_fname, 'w', format='NETCDF4')
    f.description = 'EucFACE met data, created by Martin De Kauwe'
    f.history = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date = "%s" % (datetime.datetime.now())
    f.contact = "mdekauwe@gmail.com"

    # set dimensions
    f.createDimension('time', None)
    f.createDimension('z', ndim)
    f.createDimension('y', ndim)
    f.createDimension('x', ndim)
    #f.Conventions = "CF-1.0"

    # create variables
    time = f.createVariable('time', 'f8', ('time',))
    time.units = "seconds since %s 00:00:00" % (df.Date[0])
    time.long_name = "time"
    time.calendar = "standard"

    z = f.createVariable('z', 'f8', ('z',))
    z.long_name = "z"
    z.long_name = "z dimension"

    y = f.createVariable('y', 'f8', ('y',))
    y.long_name = "y"
    y.long_name = "y dimension"

    x = f.createVariable('x', 'f8', ('x',))
    x.long_name = "x"
    x.long_name = "x dimension"

    latitude = f.createVariable('latitude', 'f8', ('y', 'x',))
    latitude.units = "degrees_north"
    latitude.missing_value = -9999.
    latitude.long_name = "Latitude"

    longitude = f.createVariable('longitude', 'f8', ('y', 'x',))
    longitude.units = "degrees_east"
    longitude.missing_value = -9999.
    longitude.long_name = "Longitude"

    SWdown = f.createVariable('SWdown', 'f8', ('time', 'y', 'x',))
    SWdown.units = "W/m^2"
    SWdown.missing_value = -9999.
    SWdown.long_name = "Surface incident shortwave radiation"
    SWdown.CF_name = "surface_downwelling_shortwave_flux_in_air"

    Tair = f.createVariable('Tair', 'f8', ('time', 'z', 'y', 'x',))
    Tair.units = "K"
    Tair.missing_value = -9999.
    Tair.long_name = "Near surface air temperature"
    Tair.CF_name = "surface_temperature"

    Rainf = f.createVariable('Rainf', 'f8', ('time', 'y', 'x',))
    Rainf.units = "mm/s"
    Rainf.missing_value = -9999.
    Rainf.long_name = "Rainfall rate"
    Rainf.CF_name = "precipitation_flux"

    Qair = f.createVariable('Qair', 'f8', ('time', 'z', 'y', 'x',))
    Qair.units = "kg/kg"
    Qair.missing_value = -9999.
    Qair.long_name = "Near surface specific humidity"
    Qair.CF_name = "surface_specific_humidity"

    Wind = f.createVariable('Wind', 'f8', ('time', 'z', 'y', 'x',))
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

    CO2 = f.createVariable('CO2air', 'f8', ('time', 'z', 'y', 'x',))
    CO2.units = "ppm"
    CO2.missing_value = -9999.
    CO2.long_name = ""
    CO2.CF_name = ""

    LAI = f.createVariable('LAI', 'f8', ('time', 'y', 'x'))
    LAI.setncatts({'long_name': u"Leaf Area Index",})

    vcmax = f.createVariable('vcmax', 'f8', ('y', 'x'))
    ejmax = f.createVariable('ejmax', 'f8', ('y', 'x'))
    g1 = f.createVariable('g1', 'f8', ('y', 'x'))

    elevation = f.createVariable('elevation', 'f8', ('y', 'x',))
    elevation.units = "m" ;
    elevation.missing_value = -9999.
    elevation.long_name = "Site elevation above sea level" ;

    reference_height = f.createVariable('reference_height', 'f8', ('y', 'x',))
    reference_height.units = "m"
    reference_height.missing_value = -9999.
    reference_height.long_name = "Crane height"

    # write data to file
    x[:] = ndim
    y[:] = ndim
    z[:] = ndim
    time[:] = times
    latitude[:] = -33.617778 # Ellsworth 2017, NCC
    longitude[:] = 150.740278 # Ellsworth 2017, NCC

    SWdown[:,0,0] = (df.PAR.values * PAR_2_SW).reshape(n_timesteps, ndim, ndim)
    Tair[:,0,0,0] = (df.TAIR.values + DEG_2_KELVIN).reshape(n_timesteps,
                                                            ndim, ndim, ndim)
    Rainf[:,0,0] = df.PPT.values.reshape(n_timesteps, ndim, ndim)
    qa_vals = convert_rh_to_qair(df.RH.values, df.TAIR.values, df.PRESS.values)
    Qair[:,0,0,0] = qa_vals.reshape(n_timesteps, ndim, ndim, ndim)
    Wind[:,0,0,0] = df.WIND.values.reshape(n_timesteps, ndim, ndim, ndim)
    PSurf[:,0,0] = df.PRESS.values.reshape(n_timesteps, ndim, ndim)
    lw = estimate_lwdown(df.TAIR.values + DEG_2_KELVIN, df.RH.values)
    LWdown[:,0,0] = lw.reshape(n_timesteps, ndim, ndim)
    if co2_conc == "amb":
        CO2[:,0,0] = df["Ca.A"].values.reshape(n_timesteps, ndim, ndim, ndim)
        vcmax[:] = 86.1425919e-6
        ejmax[:] = 138.4595736e-6
    elif co2_conc == "ele":
        CO2[:,0,0] = df["Ca.E"].values.reshape(n_timesteps, ndim, ndim, ndim)
        vcmax[:] = 81.70591263e-6
        ejmax[:] = 135.8062907e-6
    elevation[:] = 23.0 # Ellsworth 2017, NCC
    reference_height[:] = 35.0 # setting this to crane height
    LAI[:,0,0] = df.lai.values.reshape(n_timesteps, ndim, ndim)
    g1[:] = 3.8
    f.close()

def convert_rh_to_qair(rh, tair, press):
    """
    Converts relative humidity to specific humidity (kg/kg)

    Params:
    -------
    tair : float
        deg C
    press : float
        pa
    rh : float
        %
    """

    # Sat vapour pressure in Pa
    esat = calc_esat(tair)

    # Specific humidity at saturation:
    ws = 0.622 * esat / (press - esat)

    # specific humidity
    qair = (rh / 100.0) * ws

    return qair

def calc_esat(tair):
    """
    Calculates saturation vapour pressure

    Params:
    -------
    tair : float
        deg C

    Reference:
    ----------
    * Jones (1992) Plants and microclimate: A quantitative approach to
    environmental plant physiology, p110
    """

    esat = 613.75 * np.exp(17.502 * tair / (240.97 + tair))

    return esat


def estimate_lwdown(tairK, rh):
    """
    Synthesises downward longwave radiation based on Tair RH

    Reference:
    ----------
    * Abramowitz et al. (2012), Geophysical Research Letters, 39, L04808

    """
    zeroC = 273.15

    sat_vapress = 611.2 * np.exp(17.67 * ((tairK - zeroC) / (tairK - 29.65)))
    vapress = np.maximum(5.0, rh) / 100. * sat_vapress
    lw_down = 2.648 * tairK + 0.0346 * vapress - 474.0

    return lw_down


if __name__ == "__main__":

    #in_fname = "raw_data/eucdata.csv"
    in_fname = "raw_data/euc_met.csv"
    for co2_conc in ["amb", "ele"]:
        out_fname = "EucFACE_met_%s.nc" % (co2_conc)
        main(in_fname, out_fname, co2_conc=co2_conc)
