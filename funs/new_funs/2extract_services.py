import warnings
import datetime
import os
import os.path
import gc
import numpy
import math
import statistics
import rasterstats
import rasterio
import rasterio.mask
import subprocess
import json
import traceback
import shutil
from tqdm import tqdm
#import utm_zone
import numpy as np
import osmnx as ox
import networkx as nx
import pandas as pd
import geopandas as gpd
import shapely.geometry
from shapely.geometry import LineString, Point
from shapely.ops import unary_union, transform
import shapely.ops
import topojson

from  frame1.isochrones import *
from  frame1.get_service_locations import *
from  frame1.gtfs_parser import *
from  frame1.prep_bike_osm import *
from  frame1.prep_pop_ghsl import *
# import access
#import summarize_ttm

import pdb
import time


# hdc = 7277 # sao paulo
# hdc = 2007 # LA
# hdc = 8099 # ny

def extract_services(hdc)

    folder_name = 'cities_out'+'/ghsl_region_'+str(hdc)+'/'

    testing_services = []
    testing_services =  ['healthcare', 'schools', 'libraries', 'bikeshare']
            
    service_point_locations={}
    handler = ServiceHandler()
    handler.apply_file(folder_name+'temp/city.osm.pbf', locations=True)
    # extract the services
    for service in testing_services:
      # service = testing_services[1]
        coords = handler.locationlist[service]
        service_point_locations[service] = gpd.GeoDataFrame(
            geometry = [Point(coord) for coord in coords],
            crs=4326)
    citywide_carfree = handler.carfreelist
    
    
    if 'pnft' in to_test:
        testing_services.append('pnft')
        freq_stops, gtfs_wednesdays = get_frequent_stops(
        # freq_stops, gtfs_wednesdays = gtfs_parser.get_frequent_stops(
            boundaries, 
            folder_name, 
            headway_threshold)
        service_point_locations['pnft'] = freq_stops
        gtfs_filenames = [file for file in os.listdir(folder_name+'temp/gtfs/') if file[-4:] == '.zip']
        service_point_locations['pnft'].to_file(f'frame1/new_funs/data/pnft.gpkg', driver='GPKG')
    
    # export
    service_point_locations['healthcare'].to_file(f'frame1/new_funs/data/healthcare.gpkg', driver='GPKG')
    service_point_locations['schools'].to_file(f'frame1/new_funs/data/schools.gpkg', driver='GPKG')
    service_point_locations['libraries'].to_file(f'frame1/new_funs/data/libraries.gpkg', driver='GPKG')
    service_point_locations['bikeshare'].to_file(f'frame1/new_funs/data/bikeshare.gpkg', driver='GPKG')


        
    
                                
    
