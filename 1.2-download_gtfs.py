import requests
import csv
import os

import gtfs_kit as gk
import zipfile
import wget
import io

import subprocess
import datetime
import shutil
import string
import os

import pdb
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Point

from zipfile import ZipFile

from funs.get_jurisdictions import get_jurisdictions

from datetime import datetime
import numpy
import math

def get_GTFS_from_mobility_database(hdc, sources_loc='input_data/gtfs/mobilitydata_feeds_202502.csv'):
    
  # Get the current date and format it as yearmonth
  year = "2024"
  # year = datetime.now().strftime("%Y")
  today_yearmonth = "202412"
  # today_yearmonth = datetime.now().strftime("%Y%m")


  # collect boundaries
  # boundaries = get_jurisdictions(hdc)
  boundaries = gpd.read_file(f'/media/kauebraga/data/pedestriansfirst/cities_out/ghsl_region_{hdc}/debug/analysis_areas.gpkg')    
  #Let's make sure to buffer this to include peripheral roads etc for routing
  boundaries=boundaries.unary_union
  boundaries = shapely.convex_hull(boundaries)
  boundaries_latlon = gpd.GeoDataFrame(geometry = [boundaries])
  boundaries_latlon.crs = {'init':'epsg:4326'}
  longitude = round(numpy.mean(boundaries_latlon.geometry.centroid.x),10)
  utm_zone = int(math.floor((longitude + 180) / 6) + 1)
  utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
  boundaries_utm = boundaries_latlon.to_crs(utm_crs)
  boundaries_utm.geometry = boundaries_utm.geometry.buffer(1000)
  boundaries_latlon = boundaries_utm.to_crs(epsg=4326)
  boundaries_mw = boundaries_utm.to_crs("ESRI:54009")
  boundaries = boundaries_latlon.geometry.unary_union
  
  if not os.path.exists(sources_loc):
      url = 'https://bit.ly/catalogs-csv'
      r = requests.get(url, allow_redirects=True)  # to get content after redirection
      with open(sources_loc, 'wb') as f:
          f.write(r.content)
  sources = gpd.read_file(sources_loc)
  # remove deprecated?
  sources = sources[sources.status != "deprecated"]
  filenames = []
  for idx in list(sources.index):
      # idx = 2322
      if not sources.loc[idx,'bounding_box.minimum_longitude'] == '':
          sources.loc[idx,'geometry'] = shapely.geometry.box(
              float(sources.loc[idx,'bounding_box.minimum_longitude']),
              float(sources.loc[idx,'bounding_box.minimum_latitude']),
              float(sources.loc[idx,'bounding_box.maximum_longitude']),
              float(sources.loc[idx,'bounding_box.maximum_latitude']),
              )
          if sources.loc[idx,'geometry'].intersects(boundaries):
              overlap = sources.loc[idx,'geometry'].intersection(boundaries)
              if overlap.area * 1000 > sources.loc[idx,'geometry'].area:
                  url = sources.loc[idx,'urls.latest']
                  name = str(idx)+sources.loc[idx,'provider']
                  if sources.loc[idx,'name'] != '':
                      name = name+'_'+ sources.loc[idx,'name']
                  name = name.translate(str.maketrans('', '', string.punctuation))
                  name = name[:50]
                  name = name.replace(' ','_')
                  if name != '' and url != '':
                      filename=f'input_data/gtfs/{year}/gtfs_{hdc}_{today_yearmonth}_{name}.zip'
                      filenames.append(filename)
                      r = requests.get(url, allow_redirects=True)  # to get content after redirection
                      with open(filename, 'wb') as f:
                          f.write(r.content)
  return filenames


# apply function
hdcs = pd.read_csv('input_data/pbf/pbf_hdc_country.csv', dtype={'hdc' : str})
hdcs_to_test = hdcs['hdc'].tolist()

for hdc in hdcs_to_test:
    get_GTFS_from_mobility_database(hdc)
