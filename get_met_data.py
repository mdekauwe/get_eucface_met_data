#!/usr/bin/env python

"""
Get met data from the EucFACE experiment so that we can create a forcing file
for CABLE. This requires Gerry's hiev package to be installed first
(https://gdevine.github.io/hievpy/)

That's all folks.

"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (04.08.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import xarray as xr
import numpy as np

import hievpy as hp

def main():

    data = "raw_data"
    if not os.path.exists(data):
        os.makedirs(data)

    # get out code so we can access hiev
    f = open("hiev_api_token.txt", "r")
    api_token = f.read().strip()

    results = hp.search(api_token, experiments=['31'],
                        upload_from_date="2013-01-01")

    for fname in results:
        print(fname)
        hp.download(api_token, fname, path=data)


if __name__ == "__main__":

    main()
