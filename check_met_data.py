#!/usr/bin/env python
"""
Check met data looks sensible
"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import xarray as xr

__author__  = "Martin De Kauwe"
__version__ = "1.0 (22.04.2016)"
__email__   = "mdekauwe@gmail.com"

def main():

    fname = "EucFACE_met_amb.nc"
    dsa = xr.open_dataset(fname)
    fname = "EucFACE_met_ele.nc"
    dse = xr.open_dataset(fname)

    vars_to_keep = ['CO2']
    dfa = dsa[vars_to_keep].squeeze(dim=["x","y"], drop=True).to_dataframe()
    dfe = dse[vars_to_keep].squeeze(dim=["x","y"], drop=True).to_dataframe()
    co2_a = dfa.CO2.resample('D').mean()
    co2_e = dfe.CO2.resample('D').mean()

    plt.plot(co2_e)
    plt.plot(co2_a)
    #plt.plot(co2_e/co2_a)
    plt.show()

if __name__ == "__main__":

    main()
