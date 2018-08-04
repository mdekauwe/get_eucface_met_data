#!/usr/bin/env python

"""

Requires Gerry's hiev package to first be installed
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

import hievpy

def main():

    f = open("hiev_api_token.txt", "r")
    api_token = f.read().strip()
    print(api_token)


if __name__ == "__main__":

    main()
