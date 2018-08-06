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
    ds = xr.open_dataset(fname)

    x = ds.SWdown[:,0,0].values
    plt.plot(x)
    plt.show()

if __name__ == "__main__":

    main()
