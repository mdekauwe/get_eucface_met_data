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

from download_and_unzip_files import get_files

# get out code so we can access hiev
f = open("hiev_api_token.txt", "r")
api_token = f.read().strip()
url = 'https://hiev.uws.edu.au/data_files/api_search'
experiment_ids = [43]
filenames = 'BMS_S39_.*\.zip$'
get_files(api_token, url, experiment_ids, filenames)
