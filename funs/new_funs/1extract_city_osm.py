import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

from frame1.pedestriansfirst_mod import *
from frame1.get_jurisdictions import get_jurisdictions
from frame1.get_number_jurisdictions import *
from frame1.prep_poly import *
import subprocess
import os
import os.path
import json
import shutil
import geojson
import shapely
import datetime
from shapely.geometry import LineString, Polygon, shape, mapping
from shapely.ops import unary_union
import geopandas as gpd
import topojson
import pandas as pd
import numpy as np
import math
import osmnx as ox
import warnings
import urllib
import zipfile
from tqdm import tqdm
import sys
import pdb
import timeit
import csv


# hdc = 7277
# hdc = 2007 # los angeles
# hdc = 8154 # fortaleza
# hdc = 8099 # ny

def extract_city_osm(hdc, minimum_portion=0.6):
  
  #1) set up directories
    if not os.path.isdir('temp/'):
        os.mkdir('temp/')
    if not os.path.isdir('cities_out/'):
        os.mkdir('cities_out/')
    if folder_prefix:
        folder_name = 'cities_out'+'/ghsl_region_'+str(hdc)+'/'
    else:
        folder_name = str(hdc)+'/'
    if not os.path.isdir(str(folder_name)):
        os.mkdir(str(folder_name))
    if not os.path.isdir(str(folder_name)+'/debug/'):
        os.mkdir(str(folder_name)+'/debug/')
    if not os.path.isdir(str(folder_name)+'/temp/'):
        os.mkdir(str(folder_name)+'/temp/')
        
    
    analysis_areas = get_jurisdictions(
        hdc, 
        minimum_portion=minimum_portion
        )
        
    analysis_areas.to_file(f'{folder_name}/debug/analysis_areas.gpkg', driver='GPKG')
    
    #If we're going to do any access-based indicators
    #Let's make sure to buffer this to include peripheral roads etc for routing
    total_poly_latlon=analysis_areas.unary_union
    total_poly_latlon = shapely.convex_hull(total_poly_latlon)
    
    gpd.GeoDataFrame(geometry=[total_poly_latlon], crs=4326).to_file(f'{folder_name}/debug/area_for_osm_extract.gpkg', driver='GPKG')
    
    prep_from_poly(hdc, total_poly_latlon, folder_name, boundary_buffer = 2000)
