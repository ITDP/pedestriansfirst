import subprocess
import fiona
import os
import os.path
import json
import shutil
import geojson
import shapely
from shapely.geometry import LineString, Polygon, shape, mapping
from shapely.ops import unary_union
import geopandas as gpd
import topojson
import pandas as pd
import numpy
import math
import osmnx as ox
import warnings

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
    ucdb = gpd.read_file('input_data/old_ghsl/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp')
    ucdb.index =  ucdb['ID_HDC_G0']
    boundaries = ucdb.loc[hdc,'geometry']
    name = ucdb.loc[hdc,'UC_NM_MN']
    return boundaries, name

def prep_from_poly(poly, folder_name, boundary_buffer = 500):
    #save city geometry so that I can take an extract from planet.pbf within it
    #return True if overpass will be needed (planet.pbf unavailable), False otherwise
    folder_name = str(folder_name)
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
    for subfolder in ['/debug/','/temp/','/temp/access/','/geodata/']:
        subfolder = folder_name+subfolder
        if not os.path.isdir(subfolder):
            os.mkdir(subfolder)
    
    bound_latlon = gpd.GeoDataFrame(geometry = [poly], crs=4326)
    if boundary_buffer > 0:
        longitude = round(numpy.mean(bound_latlon.geometry.centroid.x),10)
        utm_zone = int(math.floor((longitude + 180) / 6) + 1)
        utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
        bound_utm = bound_latlon.to_crs(utm_crs)
        bound_utm.geometry = bound_utm.geometry.buffer(boundary_buffer)
        bound_latlon = bound_utm.to_crs(epsg=4326)
    geom_in_geojson = geojson.Feature(geometry=bound_latlon.geometry.unary_union, properties={})
    with open(folder_name+'temp/boundaries.geojson', 'w') as out:
        out.write(json.dumps(geom_in_geojson))
    #take extract from planet.pbf
    if os.path.exists('input_data/planet-latest.osm.pbf'):
        if not os.path.exists('{}/temp/city.pbf'.format(str(folder_name))):
            command = f"osmium extract input_data/planet-latest.osm.pbf -p {str(folder_name)}/temp/boundaries.geojson -s complete_ways -v -o {str(folder_name)}/temp/city.pbf"
            print(command)
            subprocess.check_call(command.split(' '))
        command = f"osmconvert {str(folder_name)}temp/city.pbf -o={str(folder_name)}temp/city.o5m"
        print(command)
        subprocess.check_call(command.split(' '))
        command = f'osmfilter {str(folder_name)}temp/city.o5m --keep="highway=" -o={str(folder_name)}temp/cityhighways.o5m'
        print(command)
        subprocess.check_call(command, shell=True)
        #todo -- read both bikeways and walkways direct from a patch'd cityhighways.osm; do walking/cycling selection logic in here.
        command = [f'osmfilter {str(folder_name)}temp/cityhighways.o5m --drop="area=yes highway=link =motor =proposed =construction =abandoned =platform =raceway service=parking_aisle =driveway =private foot=no" -o={str(folder_name)}temp/citywalk.o5m']
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
    calctime = pedestriansfirst.spatial_analysis(poly, hdc, name, folder_name, **kwargs)
    results = pedestriansfirst.calculate_indicators(poly, folder_name)
    results['calctime'] = calctime
    print(results)
    
def from_id_osm(osmid, folder_prefix = 'cities_out', boundary_buffer = 500, kwargs = {}):
    if folder_prefix:
        folder_name = folder_prefix+'/osm_'+str(osmid)
    else:
        folder_name = str(osmid)
    poly, name = poly_from_osm_cityid(osmid)
    overpass = prep_from_poly(poly, folder_name, boundary_buffer)
    kwargs['overpass'] = overpass
    calctime = pedestriansfirst.spatial_analysis(poly, osmid, name, folder_name, **kwargs)
    results = pedestriansfirst.calculate_indicators(poly, folder_name)
    results['calctime'] = calctime
    print(results)

def get_pop_ghsl(city):
    return city['properties']['P15']

#TODO: make this cut out water!
#TODO -- consider customizing this for the USA? an extra buffer or something?
def get_jurisdictions(poly_latlon,
                      minimum_portion = 0.6, #portion of a jurisdiction that has to be within the poly
                      level_min_mean_area = 5,# min size in km2 for the mean area of a unit at an admin_level
                      level_min_coverage = .2, #min coverage of an admin_level of the poly_latlon
                      buffer = 1000, #in m
                      ): 
    print('getting sub-jurisdictions')
    poly_utm_gdf = ox.projection.project_gdf(gpd.GeoDataFrame(geometry=[poly_latlon], crs=4326))
    poly_utm_gdf = poly_utm_gdf.buffer(buffer)
    poly_utm = poly_utm_gdf.unary_union
    poly_latlon = poly_utm_gdf.to_crs(4326).unary_union
    admin_lvls = [str(x) for x in range(4,11)]
    jurisdictions_latlon = ox.geometries_from_polygon(poly_latlon, tags={'admin_level':admin_lvls})
    jurisdictions_latlon = jurisdictions_latlon.loc[('relation',)]
    print(f'found {len(jurisdictions_latlon)} on first pass')
    jurisdictions_utm = ox.projection.project_gdf(jurisdictions_latlon)
    jurisdictions_clipped_utm = jurisdictions_utm.intersection(poly_utm)
    selection = (jurisdictions_clipped_utm.area / jurisdictions_utm.area) > minimum_portion
    select_jurisdictions_utm = jurisdictions_utm[selection]
    print(f'found {len(select_jurisdictions_utm)} with {minimum_portion} inside area')
    select_jurisdictions_latlon = select_jurisdictions_utm.to_crs(4326)
    if len(select_jurisdictions_latlon) > 0:
        total_boundaries_latlon = unary_union([select_jurisdictions_latlon.unary_union, poly_latlon])
        total_boundaries_utm = unary_union([select_jurisdictions_utm.unary_union, poly_utm])
    else:
        total_boundaries_latlon = poly_latlon
        total_boundaries_utm = poly_utm
        
    #now get all jurisdictions within total_boundaries
    jurisdictions_latlon = ox.geometries_from_polygon(total_boundaries_latlon, tags={'admin_level':admin_lvls})
    jurisdictions_latlon = jurisdictions_latlon.loc[('relation',)]
    print(f'found {len(jurisdictions_latlon)} on second pass')
    jurisdictions_utm = ox.projection.project_gdf(jurisdictions_latlon)
    jurisdictions_clipped_utm = jurisdictions_utm.intersection(total_boundaries_utm)
    selection = (jurisdictions_clipped_utm.area / jurisdictions_utm.area) > 0.95
    select_jurisdictions_utm = jurisdictions_utm[selection]
    print(f'found {len(select_jurisdictions_utm)} with 0.95 inside total area')
    selected_levels = []
    for admin_level in select_jurisdictions_utm.admin_level.unique():
        selection = select_jurisdictions_utm[select_jurisdictions_utm.admin_level==admin_level]
        if selection.area.mean() >= level_min_mean_area*1000000:
            if selection.unary_union.area >= (level_min_coverage * poly_utm.area):
                selected_levels.append(admin_level)
    final_jurisdictions_utm = select_jurisdictions_utm[select_jurisdictions_utm.admin_level.isin(selected_levels)]
    final_jurisdictions_latlon = final_jurisdictions_utm.to_crs(4326)
    print(f'found {len(final_jurisdictions_latlon)} in acceptable admin levels {selected_levels}')

    return final_jurisdictions_latlon

def regional_analysis(hdc, 
                      folder_prefix = 'cities_out', 
                      minimum_portion=0.6,
                      analyze=True,
                      summarize=True,
                      simplification=0.001 #toposimplification factor
                      ):
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
    ghsl_boundaries, name = poly_from_ghsl_hdc(hdc)
    jurisdictions_latlon = get_jurisdictions(
        ghsl_boundaries, 
        minimum_portion=minimum_portion
        )
    total_poly_latlon=jurisdictions_latlon.unary_union
    
    analysis_areas = gpd.GeoDataFrame()
    new_id = 0
    analysis_areas.loc[new_id,'name'] = name
    analysis_areas.loc[new_id, 'geometry'] = ghsl_boundaries
    analysis_areas.loc[new_id, 'hdc'] = hdc
    analysis_areas.crs=4326
    new_id += 1
    for osmid in jurisdictions_latlon.index:
        analysis_areas.loc[new_id,'osmid'] = osmid
        analysis_areas.loc[:,'geometry'].loc[new_id] = jurisdictions_latlon.loc[osmid,'geometry']
        #the above hack is necessary because sometimes geometry is a multipolygon
        for attr in ['name','admin_level']:
            analysis_areas.loc[new_id,attr] = jurisdictions_latlon.loc[osmid,attr]
        new_id += 1
    
    analysis_areas.to_file(f'{folder_name}/debug/analysis_areas.gpkg', driver='GPKG')
    
    prep_from_poly(total_poly_latlon, folder_name)
    
    #now actually call the functions and save the results
    if analyze == True:
        calctime = pedestriansfirst.spatial_analysis(
            total_poly_latlon,
            hdc,
            name,
            folder_name=folder_name,
            )
    else:
        calctime = 0
    
    if summarize == True:
        all_results = pd.DataFrame()
        print('getting Urban Centre results')
        agglomeration_results = pedestriansfirst.calculate_indicators(
            ghsl_boundaries,
            folder_name=folder_name)
        for result in agglomeration_results.keys():
            all_results.loc['agglomeration', result] = agglomeration_results[result] 
        all_results.loc['agglomeration', 'geospatial_calctime'] = calctime
        for area_id in analysis_areas.index:
            print(f'getting results for {analysis_areas.loc[area_id,"name"]}')
            juri_results = pedestriansfirst.calculate_indicators(
                analysis_areas.loc[area_id,'geometry'],
                folder_name=folder_name)
            all_results.loc[area_id, 'name'] = analysis_areas.loc[area_id,'name']
            all_results.loc[area_id, 'osmid'] = analysis_areas.loc[area_id,'osmid']
            all_results.loc[area_id, 'hdc'] = analysis_areas.loc[area_id,'hdc']
            for result in juri_results.keys():
                analysis_areas.loc[area_id, result] = juri_results[result] 
                all_results.loc[area_id, result] = juri_results[result] 
            all_results.loc[area_id, 'geospatial_calctime'] = calctime
        
        #spatially simplify while preserving topology
        topo = topojson.Topology(analysis_areas, prequantize=True)
        analysis_areas = topo.toposimplify(0.001).to_gdf()
        
        all_results.to_csv(f'{folder_name}indicator_values.csv')
        analysis_areas.to_file(f'{folder_name}indicator_values.gpkg',driver='GPKG')
    
    #import pdb; pdb.set_trace()
        

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
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    ox.utils.config(log_console = False)
    
    #test PNRT in DC
    regional_analysis(12080)
    
    
    test_cities = [
        #154,
        #621,
        #5134,
        #8050,
        #855,
        #14,
        1105,
        634,
        3902,
        1406,
        1709,
        88,
        3562,
        11862,
        4427,
        4172,
        4608,
        10076,
        350,
        200,
        9691,
        1372,
        7041,
        6522,
        1445,
        1361,
        931,
        13039,
        ]
    
    for cityid in test_cities:
        regional_analysis(cityid)

    pass


