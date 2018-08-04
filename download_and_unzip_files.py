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
import json
import urllib.request
from urllib.request import urlopen
from urllib.request import Request
import zipfile
from datetime import datetime
import hievpy as hp

def get_files(api_token, url, experiment_ids, filenames, download_dir='data'):

    request_data = json.dumps({'auth_token': api_token,
                               'experiments':experiment_ids,
                               'filename':filenames}).encode("utf-8")
    request_headers = {'Content-Type' : 'application/json; charset=UTF-8',
                       'X-Accept': 'application/json'}
    req = Request(url, request_data, request_headers)
    response = urlopen(req)
    js = json.load(response)

    print( 'Number of files = %s\n' % (len(js)) )

    if len(js):
        where = os.path.join(os.path.dirname(__file__))
        dest_dir = os.path.join(where, download_dir)
        unzipped_dir = os.path.join(os.path.dirname(__file__),
                                    download_dir, 'unzipped')
        if not os.path.exists(unzipped_dir):
            os.makedirs(unzipped_dir)

        for entry in js:
            if not os.path.isfile(os.path.join(dest_dir, entry['filename'])):

                download_url = entry['url']+'?'+'auth_token=%s' %api_token
                request = Request(download_url)
                f = urlopen(request)
                path = os.path.join(dest_dir, entry['filename'])
                with open(path, 'wb') as local_file:
                    local_file.write(f.read())
                local_file.close()

                # If zipped, unzip it
                if entry['filename'].endswith(".zip"):
                    zip_fname = os.path.splitext(entry['filename'])[0]
                    if not os.path.join(unzipped_dir, zip_fname):
                        os.makedirs(unzipped_dir, zip_fname)

                    path = os.path.join(dest_dir, entry['filename'])
                    zfile = zipfile.ZipFile(path)
                    zfile.extractall(os.path.join(unzipped_dir, zipname))
    else:
        raise('No matching files!')
