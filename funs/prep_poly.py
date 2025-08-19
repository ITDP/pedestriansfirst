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
import re


def poly_from_osm_cityid(osmid):
    #return shapely polygon
    admin_area = ox.geocode_to_gdf("R"+str(osmid), by_osmid=True)
    boundaries = admin_area.geometry[0]
    name = admin_area.display_name[0]
    return boundaries, name

# poly = total_poly_latlon


# hdc = 134
def prep_from_poly(hdc, poly, folder_name, boundary_buffer = 500, current_year = 2024):
    #save city geometry so that I can take an extract from planet.pbf within it
    #return True if overpass will be needed (planet.pbf unavailable), False otherwise
    # open file with the equivalence of hdc and country osm file to open
    pbf_hdc_country = pd.read_csv('input_data/pbf/pbf_hdc_country.csv', dtype={'hdc' : str})
    # pbf_hdc_country = pd.read_csv('input_data/pbf/pbf_hdc_country.csv')
    # filter the hdc in question
    pbf_to_open = pbf_hdc_country[pbf_hdc_country['hdc'] == hdc]
    country_to_open = pbf_to_open[['country1']].values[0]
    country_to_open = ''.join(country_to_open)

    directory = f"input_data/pbf/{current_year}"    

    file_dir = [f for f in os.listdir(directory) if re.search(f"^{country_to_open}_{current_year}", f)]
    file_dir = [os.path.join(directory, f) for f in file_dir]
    file_dir = ''.join(str(x) for x in file_dir)
    
    # file_dir = f"input_data/pbf/{current_year}/{country_to_open}_{current_year}.osm.pbf"
    
    # add year
    # re.sub(r"(pbf/)", r"\1" + f"{current_year}/" + "/", file_dir)
    
    folder_name = str(folder_name)
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
    for subfolder in ['/debug/','/temp/','/temp/access/','/geodata/','/geodata/access/']:
        subfolder = folder_name+subfolder
        if not os.path.isdir(subfolder):
            os.mkdir(subfolder)
    
    bound_latlon = gpd.GeoDataFrame(geometry = [poly], crs=4326)
    if boundary_buffer > 0:
        longitude = round(np.mean(bound_latlon.geometry.centroid.x),10)
        utm_zone = int(math.floor((longitude + 180) / 6) + 1)
        utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
        bound_utm = bound_latlon.to_crs(utm_crs)
        bound_utm.geometry = bound_utm.geometry.buffer(boundary_buffer)
        bound_latlon = bound_utm.to_crs(epsg=4326)
    geom_in_geojson = geojson.Feature(geometry=bound_latlon.geometry.unary_union, properties={})
    with open(folder_name+'temp/boundaries.geojson', 'w') as out:
        out.write(json.dumps(geom_in_geojson))
    #take extract from planet.pbf
    if not os.path.exists('{}/temp/city.pbf'.format(str(folder_name))):
        command = f"osmium extract {file_dir} -p {str(folder_name)}/temp/boundaries.geojson -s complete_ways -v -o {str(folder_name)}/temp/city_{current_year}.osm.pbf --overwrite"
        print(command)
        subprocess.check_call(command.split(' '))
    if not os.path.exists(f"{str(folder_name)}temp/cityhighways_{current_year}.o5m"):
        # command = f"osmosis --read-pbf {str(folder_name)}temp/city.pbf --write-xml {str(folder_name)}temp/city.osm"
        # # command = f"osmconvert {str(folder_name)}temp/city.pbf -o={str(folder_name)}temp/city.o5m"
        # print(command)
        # subprocess.check_call(command.split(' '))
        command = f"osmium tags-filter {str(folder_name)}/temp/city_{current_year}.osm.pbf w/highway -o {str(folder_name)}temp/cityhighways_{current_year}.osm --overwrite",
        # command = f"osmium tags-filter {str(folder_name)}/temp/city.osm.pbf nwr highway,type,cycleway,bicycle,cycleway:left,cycleway:right,cycleway:both,area,service,foot,bridge,tunnel,oneway,lanes,ref,name,maxspeed,access,landuse,width,est_width,junction -o {str(folder_name)}temp/cityhighways.osm --overwrite",
        # command = f'osmfilter {str(folder_name)}temp/city.osm --keep="highway=" --keep-tags="all type= highway= cycleway= bicycle= cycleway:left= cycleway:right= cycleway:both= area= service= foot= bridge= tunnel= oneway= lanes= ref= name= maxspeed= access= landuse= width= est_width= junction=" -o={str(folder_name)}temp/cityhighways.osm'
        # command = f'osmfilter {str(folder_name)}temp/city.o5m --keep="highway=" --keep-tags="all type= highway= cycleway= bicycle= cycleway:left= cycleway:right= cycleway:both= area= service= foot= bridge= tunnel= oneway= lanes= ref= name= maxspeed= access= landuse= width= est_width= junction=" -o={str(folder_name)}temp/cityhighways.o5m'
        print(command)
        subprocess.check_call(command, shell=True)
        # do only for bikes
        command = f"osmium tags-filter {str(folder_name)}/temp/city_{current_year}.osm.pbf w/highway -o {str(folder_name)}temp/cityhighways_{current_year}.osm --overwrite",
        subprocess.check_call(command, shell=True)
        #todo -- read both bikeways and walkways direct from a patch'd cityhighways.osm; do walking/cycling selection logic in here.
        command = [f'osmfilter {str(folder_name)}temp/cityhighways_{current_year}.osm --drop="area=yes highway=link =motor =proposed =construction =abandoned =platform =raceway service=parking_aisle =driveway =private foot=no" -o={str(folder_name)}temp/citywalk_{current_year}.osm']
        print(command)
        subprocess.check_call(command, shell=True)
    return False
