import subprocess
import fiona
import os
import os.path
import json
import shutil
import geojson
import shapely
from shapely.geometry import LineString, Polygon, shape, mapping
import geopandas as gpd
import pandas as pd
import numpy
import math
import osmnx as ox

import pdb

import pedestriansfirst

def poly_from_osm_cityid(osmid):
    #return shapely polygon
    admin_area = ox.geocode_to_gdf("R"+str(osmid), by_osmid=True)
    boundaries = admin_area.geometry[0]
    name = admin_area.display_name[0]
    return boundaries, name

def poly_from_ghsl_hdc(hdc):
    #select city from ID number, return shapely polygon
    with fiona.open('input_data/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp','r') as ucdb:
        for city in ucdb:
            if int(city['properties']['ID_HDC_G0']) == int(hdc):
                target = city
    boundaries = shape(target['geometry'])
    name = target['properties']['UC_NM_MN']
    return boundaries, name

def prep_from_poly(poly, folder_name, boundary_buffer = 500):
    #save city geometry so that I can take an extract from planet.pbf within it
    #return True if overpass will be needed (planet.pbf unavailable), False otherwise
    if not os.path.isdir(str(folder_name)):
        os.mkdir(str(folder_name))
    if not os.path.isdir(str(folder_name)+'/debug/'):
        os.mkdir(str(folder_name)+'/debug/')
    bound_latlon = gpd.GeoDataFrame(geometry = [poly], crs=4326)
    if boundary_buffer > 0:
        longitude = round(numpy.mean(bound_latlon.geometry.centroid.x),10)
        utm_zone = int(math.floor((longitude + 180) / 6) + 1)
        utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
        bound_utm = bound_latlon.to_crs(utm_crs)
        bound_utm.geometry = bound_utm.geometry.buffer(boundary_buffer)
        bound_latlon = bound_utm.to_crs(epsg=4326)
    geom_in_geojson = geojson.Feature(geometry=bound_latlon.geometry.unary_union, properties={})
    with open(folder_name+'/boundaries.geojson', 'w') as out:
        out.write(json.dumps(geom_in_geojson))
    #take extract from planet.pbf
    if os.path.exists('input_data/planet-latest.osm.pbf'):
        if not os.path.exists('{}/city.pbf'.format(str(folder_name))):
            command = f"osmium extract input_data/planet-latest.osm.pbf -p {str(folder_name)}/boundaries.geojson -s complete_ways -v -o {str(folder_name)}/city.pbf"
            print(command)
            subprocess.check_call(command.split(' '))
        command = f"osmconvert {str(folder_name)}/city.pbf -o={str(folder_name)}/city.o5m"
        print(command)
        subprocess.check_call(command.split(' '))
        command = f'osmfilter {str(folder_name)}/city.o5m --keep="highway=" -o={str(folder_name)}/cityhighways.o5m'
        print(command)
        subprocess.check_call(command, shell=True)
        #todo -- read both bikeways and walkways direct from a patch'd cityhighways.osm; do walking/cycling selection logic in here.
        command = [f'osmfilter {str(folder_name)}/cityhighways.o5m --drop="area=yes highway=link =motor =proposed =construction =abandoned =platform =raceway service=parking_aisle =driveway =private foot=no" -o={str(folder_name)}/citywalk.o5m']
        print(command)
        subprocess.check_call(command, shell=True)
        return False
    else:
        return True
    
    
def from_id_hdc(hdc, folder_prefix = 'cities_out', boundary_buffer = 500, kwargs = {}):
    if folder_prefix:
        folder_name = folder_prefix+'/ghsl_'+str(hdc)
    else:
        folder_name = str(hdc)
    poly, name = poly_from_ghsl_hdc(hdc)
    overpass = prep_from_poly(poly, folder_name, boundary_buffer)
    kwargs['overpass'] = overpass
    return pedestriansfirst.pedestrians_first(poly, hdc, name, folder_name, **kwargs)
    
def from_id_osm(osmid, folder_prefix = 'cities_out', boundary_buffer = 500, kwargs = {}):
    if folder_prefix:
        folder_name = folder_prefix+'/osm_'+str(osmid)
    else:
        folder_name = str(osmid)
    poly, name = poly_from_osm_cityid(osmid)
    overpass = prep_from_poly(poly, folder_name, boundary_buffer)
    kwargs['overpass'] = overpass
    return pedestriansfirst.pedestrians_first(poly, osmid, name, folder_name, **kwargs)

def get_pop_ghsl(city):
    return city['properties']['P15']

#all cities in descending order
def all_cities():
    with fiona.open('GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp','r') as ucdb:
        cities = list(ucdb)
    cities.sort(key=get_pop_ghsl, reverse = True)
    for city in cities:
        if os.path.exists('all_results.json'):
            with open('all_results.json','r') as in_file:
                all_results = json.load(in_file)
        else:
            all_results = {}
        if not str(city['properties']['ID_HDC_G0']) in all_results.keys():
            if not str(city['properties']['ID_HDC_G0']) == '4541': #there's one city in south sudan that doesn't work right.
                results = from_id_hdc(city['properties']['ID_HDC_G0'])
                all_results.update({city['properties']['ID_HDC_G0']:results})
                with open('all_results.json','w') as out_file:
                    json.dump(all_results, out_file)



def pnb_run():
    osmids = list(pd.read_csv('cities_for_pnb.csv')['OSM ID'])

    for osmid in osmids:
        osmid=str(int(osmid))
        print('pnb run',osmid)
        if os.path.exists('pnb_results.csv'):
            results = pd.read_csv('pnb_results.csv', index_col=0)
        else:
            pnb_results = pd.DataFrame()
        if not str(osmid) in pnb_results.columns:
            results = from_id_osm(osmid, kwargs={'to_test':['pnpb','pnab','density'], 'debug':True})
            pnb_results[osmid]=results
            pnb_results.to_csv('pnb_results.csv')
                
if __name__ == '__main__':
    results_bigmax = from_id_osm(2220789, kwargs={'to_test':['pnpb','pnab','density','h+s',], 'max_edge_length':10000, 'debug':True})
    results_lilmax = from_id_osm(2220789, kwargs={'to_test':['pnpb','pnab','density','h+s',], 'max_edge_length':100, 'debug':True})
    print('big',results_bigmax)
    print('lil',results_lilmax)
    #results_kampala = from_id_hdc(4427, kwargs={'to_test':['pnpb','pnab','density','h+s','special'], 'debug':True})
    #print('big',results_bigmax)
    #print('lil',results_lilmax)
    #print(results_kampala)

#hdcs = { #test
#'Mexico City': 154,
#        }



#for city in hdcs.keys():
#    if not os.path.exists(str(hdcs[city])+'/results.json'):
#        from_id_hdc(hdcs[city])
#    else:
#        for file in ['city.o5m','cityhighways.o5m','citywalk.o5m']:
#            if os.path.exists(str(hdcs[city])+'/'+file):
#                os.remove(str(hdcs[city])+'/'+file)

