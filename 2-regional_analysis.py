import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

from funs.pedestriansfirst import *
from funs.get_jurisdictions import get_jurisdictions
from funs.get_number_jurisdictions import *
from funs.prep_poly import *
# from funs.process_bike_new import process_bike
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



import pdb

import funs.pedestriansfirst
import re
import funs.calculate_country_ind
from funs.calculate_country_ind import *

    
useful_tags = ox.settings.useful_tags_way + ['cycleway', 'cycleway:left', 'cycleway:right', 'cycleway:both', 'bicycle']
ox.settings.log_console = True
ox.settings.use_cache = False
ox.settings.useful_tags_way = useful_tags
ox.settings.overpass_rate_limit = False

# hdc = 11480
# hdc = "02007" # los angeles
# hdc = "08154" # fortaleza
# hdc = "07277" # sao paulo
# hdc = 40
# hdc = 8099 # new york
# hdc = "01156" # trujillo
# hdc = "05472" # jakarta
# hdc = "05816" # london
# hdc = "00992" # amman
# hdc = "00213" # kohima
# hdc = "02006" # cdmx
# hdc = "08238" # Salvador
# hdc = "01590" # Amsterdam


# folder_prefix = '/media/kauebraga/data/pedestriansfirst/cities_out'; current_year=2024; minimum_portion=0.5;


def regional_analysis(hdc, to_test, folder_prefix = '/media/kauebraga/data/pedestriansfirst/cities_out', minimum_portion=0.5, prep = True, analyse=True,summarize=True,jurisdictions=True,simplification=0.001, current_year=2024, cleanup = False):
    """Run analysis scripts for a given region
    
    Uses lots (too many?) hardcoded defaults
    
    1) sets up directory structure for an agglomeration
    
    2) calls get_jurisdictions to identify analysis areas within agglomeration
    
    3) prepares OSM data
    
    4)'Analysis' creates all the geodata needed to calculate
    Atlas indicators for a given analysis area (eg., isochrone polys). 
    
    5)'Summary' actually calculates those indicators for each analysis area
    within the agglomeration
    
    
    """
    
    
    useful_tags = ox.settings.useful_tags_way + ['cycleway', 'cycleway:left', 'cycleway:right', 'cycleway:both', 'bicycle']
    ox.settings.log_console = True
    ox.settings.use_cache = False
    ox.settings.useful_tags_way = useful_tags
    ox.settings.overpass_rate_limit = False
    
    df = pd.DataFrame()
    
    # Define the CSV file path
    csv_file = f"logs/running_{hdc}.csv"
    df.to_csv(csv_file, mode='a', header=True, index=False)
    
    #1) set up directories
    if not os.path.isdir('temp/'):
        os.mkdir('temp/')
    if not os.path.isdir('cities_out/'):
        os.mkdir('cities_out/')
    if folder_prefix:
        folder_name = folder_prefix+'/ghsl_region_'+str(hdc)+'/'
    else:
        folder_name = str(hdc)+'/'
    if not os.path.isdir(str(folder_name)):
        os.mkdir(str(folder_name))
    if not os.path.isdir(str(folder_name)+'/debug/'):
        os.mkdir(str(folder_name)+'/debug/')
    if not os.path.isdir(str(folder_name)+'/temp/'):
        os.mkdir(str(folder_name)+'/temp/')
    
    if os.path.isfile(f"{folder_name}geodata/blocks/blocks_latlon_{current_year}.geojson"):
        analyze = False
        
    #2) get analysis areas
    
    
    
    if jurisdictions == True:
      analysis_areas = get_jurisdictions(
          hdc, 
          minimum_portion=minimum_portion
          )
      
      analysis_areas.to_file(f'{folder_name}/debug/analysis_areas.gpkg', driver='GPKG')
    else:
      analysis_areas = gpd.read_file(f'{folder_name}/debug/analysis_areas.gpkg')    
      

    
    
    e#If we're going to do any access-based indicators
    #Let's make sure to buffer this to include peripheral roads etc for routing
    total_poly_latlon=analysis_areas.unary_union
    total_poly_latlon = shapely.convex_hull(total_poly_latlon)   
    # total_poly_latlon.to_file(f'{folder_name}/debug/boundary.gpkg', driver='GPKG')
    gpd.GeoDataFrame(geometry=[total_poly_latlon], crs=4326).to_file(f'{folder_name}/debug/area_for_osm_extract.gpkg', driver='GPKG')
    
    #3) prepare OSM data files for the agglomeration
    if prep == True:
      start = timeit.default_timer()
      prep_from_poly(hdc, total_poly_latlon, folder_name, boundary_buffer = 2000, current_year=current_year)
      stop = timeit.default_timer()
      print('Time prep: ', stop - start)  
    else:
      print('Didnt prep')  
    
    #4) now actually call the functions and save the results
    if analyse:
      print("why?")
      start = timeit.default_timer()
      geospatial_calctime = spatial_analysis(
          boundaries = total_poly_latlon,
          id_code = hdc,
          name = analysis_areas.loc[0,'name'],
          folder_name=folder_name,
          current_year = current_year,
          to_test = to_test
          )
      stop = timeit.default_timer()
      print('Time analyze: ', stop - start)  
    else:
      print('Didnt analyze')  
      geospatial_calctime = 0
    
    #5) calculate indicator measurement for each analysis area
    if summarize == True:
        start = datetime.datetime.now()
        calculate_indicators(
            analysis_areas = analysis_areas, 
            folder_name=folder_name,
            current_year = current_year
            )
        end = datetime.datetime.now()
        analysis_areas.loc[0, 'geospatial_calctime'] = str(geospatial_calctime)
        analysis_areas.loc[0, 'summary_calctime'] = str(end - start)
        
        topo = topojson.Topology(analysis_areas, prequantize=True)
        analysis_areas = topo.toposimplify(0.001).to_gdf()
        analysis_areas.to_file(f'{folder_name}indicator_values_{current_year}.gpkg',driver='GPKG')
        
        nongeospatial_results = analysis_areas.drop('geometry', axis=1, inplace=False)
        nongeospatial_results.to_csv(f'{folder_name}indicator_values_{current_year}.csv')
        
    # Define the CSV file path
    os.remove(csv_file)
    csv_file = f"logs/finished_{hdc}.csv"
    df.to_csv(csv_file, mode='a', header=True, index=False)
    
    #clean up big files
    if cleanup:
        for cleanup_filename in ['city.o5m', 'city.pbf','cityhighways.o5m','citywalk.o5m','access/city_ltstagged.pbf']:
            if os.path.exists(f'{folder_name}/temp/{cleanup_filename}'):
                os.remove(f'{folder_name}/temp/{cleanup_filename}')
        for file in os.listdir('cache/'):
            os.remove(f'cache/{file}')
    #import pdb; pdb.set_trace()


