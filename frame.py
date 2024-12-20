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

import pedestriansfirst

def poly_from_osm_cityid(osmid):
    #return shapely polygon
    admin_area = ox.geocode_to_gdf("R"+str(osmid), by_osmid=True)
    boundaries = admin_area.geometry[0]
    name = admin_area.display_name[0]
    return boundaries, name

def prep_from_poly(poly, folder_name, boundary_buffer = 500):
    #save city geometry so that I can take an extract from planet.pbf within it
    #return True if overpass will be needed (planet.pbf unavailable), False otherwise
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
    if os.path.exists('input_data/planet-latest.osm.pbf'):
        if not os.path.exists('{}/temp/city.pbf'.format(str(folder_name))):
            command = f"osmium extract input_data/planet-latest.osm.pbf -p {str(folder_name)}/temp/boundaries.geojson -s complete_ways -v -o {str(folder_name)}/temp/city.pbf"
            print(command)
            subprocess.check_call(command.split(' '))
        if not os.path.exists(f"{str(folder_name)}temp/cityhighways.o5m"):
            command = f"osmconvert {str(folder_name)}temp/city.pbf -o={str(folder_name)}temp/city.o5m"
            print(command)
            subprocess.check_call(command.split(' '))
            command = f'osmfilter {str(folder_name)}temp/city.o5m --keep="highway=" --keep-tags="all type= highway= cycleway= bicycle= cycleway:left= cycleway:right= cycleway:both= area= service= foot= bridge= tunnel= oneway= lanes= ref= name= maxspeed= access= landuse= width= est_width= junction=" -o={str(folder_name)}temp/cityhighways.o5m'
            print(command)
            subprocess.check_call(command, shell=True)
            #todo -- read both bikeways and walkways direct from a patch'd cityhighways.osm; do walking/cycling selection logic in here.
            command = [f'osmfilter {str(folder_name)}temp/cityhighways.o5m --drop="area=yes highway=link =motor =proposed =construction =abandoned =platform =raceway service=parking_aisle =driveway =private foot=no" -o={str(folder_name)}temp/citywalk.o5m']
            print(command)
            subprocess.check_call(command, shell=True)
        return False
    else:
        return True
    
    
    
def download_ghsl(proj='mw', resolution='1000'):
    """Download population data from the Global Human Settlement Layer
    
    Can get data in mollewiede or WGS84 projections
    Hard-coded to use 2023 release
    """
    
    proj_code = {'mw':54009,'ll':4326}[proj]
    if not os.path.exists('input_data/'):
        os.mkdir('input_data/')
    if not os.path.exists(f'input_data/ghsl_data_{resolution}m_{proj}/'):
        os.mkdir('input_data/ghsl_data_{resolution}m_mw/')
    for year in tqdm(range(1975, 2031, 5)): 
        letter='E'
        if proj == 'mw':
            name = f'GHS_POP_{letter}{year}_GLOBE_R2023A_{proj_code}_{resolution}'
        if not os.path.exists(f'input_data/ghsl_data_{resolution}m_{proj}/{name}/{name}_V1_0.tif'):
            zippath, _ = urllib.request.urlretrieve(f'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/{name}/V1-0/{name}_V1_0.zip')
            with zipfile.ZipFile(zippath, "r") as f:
                f.extractall(f'input_data/ghsl_data_{resolution}m_{proj}/{name}')
    
    
#TODO: make this cut out water!
#TODO -- consider customizing this for the USA? an extra buffer or something?
def get_jurisdictions(hdc,
                      minimum_portion = 0.5, #portion of a jurisdiction that has to be within the poly
                      level_min_mean_area = 2,# min size in km2 for the mean area of a unit at an admin_level
                      level_min_coverage = .0000002, #min coverage of an admin_level of the poly_latlon
                      buffer = 4000, #in m
                      ): 
    """Get OSM data for administratrive jursidictions within an agglomeration
    """
    
    #Special call to get ITDP "Cycling Cities" where the jurisdiction 
    #is an extra distance from the agglomeration
    if hdc in [8265, 102]: #kohima, zapopan
        buffer = 40000
    
    #load Urban Centre Database, assign full names
    ucdb = gpd.read_file('input_data/ghsl/SMOD_V1s6_opr_P2023_v1_2020_labelUC_DB_release.gpkg')
    ucdb.index =  ucdb['ID_UC_G0']
    ghsl_boundaries_mw = ucdb.loc[hdc,'geometry']
    name_full = ucdb.loc[hdc,'NAME_LIST']
    all_names = name_full.split('; ')
    name = "The " + " / ".join(all_names[:3]) + ' area'
    if len(name) >= 50:
        name = "The " + " / ".join(all_names[:2]) + ' area'
        if len(name) >= 50:
            name = "The " + all_names[0] + ' area'
    name_short = "The " + all_names[0] + ' area'
    
    #convert to UTM and buffer
    poly_mw_gdf = gpd.GeoDataFrame(geometry=[ghsl_boundaries_mw], crs="ESRI:54009")
    poly_ll_gdf = poly_mw_gdf.to_crs(4326)
    ghsl_boundaries = poly_ll_gdf.unary_union
    poly_utm_gdf = ox.project_gdf(poly_ll_gdf)
    buffered_poly_utm_gdf = poly_utm_gdf.buffer(buffer)
    buffered_poly_latlon_gdf = buffered_poly_utm_gdf.to_crs(4326)
    buffered_poly_utm = buffered_poly_utm_gdf.unary_union
    buffered_poly_latlon = buffered_poly_latlon_gdf.unary_union
    
    #set up dataframe for output
    analysis_areas = gpd.GeoDataFrame()
    new_id = 0
    analysis_areas.loc[new_id, 'name'] = name
    analysis_areas.loc[new_id, 'name_short'] = name_short
    analysis_areas.loc[new_id, 'geometry'] = ghsl_boundaries
    analysis_areas.loc[new_id, 'hdc'] = hdc
    analysis_areas.loc[new_id, 'osmid'] = None
    analysis_areas.loc[new_id, 'level_name'] = 'Agglomeration'
    analysis_areas.crs=4326
    
    #there's probably a better way of doing this. I increment this 
    #counter in order to set a unique new ID for every new analysis_area
    new_id += 1
    
    #now figure out what country it's in
    #and while we're at it, set up country-specific analysis areas
    
    #get CGAZ data here, also use it for clipping coastline
    country_bounds = gpd.read_file('input_data/CGAZ/geoBoundaries_ITDPv5.gpkg')
    country_bounds.crs=4326
    country_bounds.geometry = country_bounds.make_valid()
    earth_utm = country_bounds.to_crs(crs = poly_utm_gdf.crs)
    #get land within 10km
    area_for_land_ll = buffered_poly_utm_gdf.buffer(100000).to_crs(4326).unary_union
    
    try:
        country_land_sum = country_bounds.intersection(area_for_land_ll).area.sum()
    except: 
        pdb.set_trace() #debug
    if country_land_sum >= area_for_land_ll.area * 0.95:
        nearby_land_gdf_utm = buffered_poly_utm_gdf.buffer(100000)
        nearby_land_gdf_ll = buffered_poly_utm_gdf.to_crs(4326)
    else:
        nearby_land_gdf_ll = gpd.clip(country_bounds, area_for_land_ll)
        nearby_land_gdf_utm = nearby_land_gdf_ll.to_crs(buffered_poly_utm_gdf.crs)
            
    country_overlaps = country_bounds.overlay(poly_ll_gdf, how='intersection')
    country_overlaps['land_area'] = country_overlaps.to_crs(poly_utm_gdf.crs).area
    main_country = country_overlaps.sort_values('land_area', ascending=False).iloc[0]['shapeGroup']
    
    # we create analysis areas for "the portion of an agglomeration within each 
    # country
    # May no longer be needed with new UCDB data (coming Nov 2024)
    countries = list(country_overlaps.sort_values('shapeGroup', ascending=True).shapeGroup)
    agglomeration_country_name = ' / '.join(countries)
    analysis_areas.loc[0, 'agglomeration_country_name'] = agglomeration_country_name

    # add country-specific analysis areas
    for idx in country_overlaps.index:
        analysis_areas.loc[new_id,'country'] = country_overlaps.loc[idx,'shapeGroup']
        analysis_areas.loc[:,'geometry'].loc[new_id] = country_overlaps.loc[idx,'geometry']
        new_id += 1
    
    #first, for Brazil only, we check which 'metropolitan areas' it's in
    #based on data from ITDP Brazil
    #1. get the area that includes the centroid of the agglomeration
    #2. get any other areas that are at least ?30%? covered by the agglomeration
    brazil_metros = gpd.read_file('input_data/country_specific_zones/brazil_selected_metro_areas.gpkg')
    brazil_metros.to_crs(4326)
    if main_country == 'BRA':
        print("We're in Brazil :)")
        brazil_metros_utm = brazil_metros.to_crs(buffered_poly_utm_gdf.crs)
        overlap = brazil_metros_utm.intersection(buffered_poly_utm)
        selection = (overlap.area / brazil_metros_utm.area) > 0.3
        if len(brazil_metros[selection]) == 0:
            selection = brazil_metros_utm.intersects(buffered_poly_utm.centroid)
        select_metro_areas_utm = brazil_metros_utm[selection]
        select_metro_areas_latlon = select_metro_areas_utm.to_crs(4326)
        for area in select_metro_areas_latlon.iterrows():
            analysis_areas.loc[new_id, 'name_long'] = area[1].name_muni
            analysis_areas.loc[new_id, 'geometry'] = area[1].geometry
            analysis_areas.loc[new_id, 'hdc'] = None
            analysis_areas.loc[new_id, 'osmid'] = None
            analysis_areas.loc[new_id, 'level_name'] = 'Brazilian Metro Areas'
            new_id += 1
        #and union it with the poly_latlons, mostly so we get jurisdictions 
        #inside brasilia
        buffered_poly_latlon = unary_union([
            select_metro_areas_latlon.unary_union,
            buffered_poly_latlon
            ])
        buffered_poly_utm = unary_union([
            select_metro_areas_utm.unary_union,
            buffered_poly_utm
            ])
        buffered_poly_utm_gdf = gpd.GeoDataFrame(geometry = [buffered_poly_utm], crs = poly_utm_gdf.crs)
    
        
    #STEP 1, get all the sub-jusisdictions at least minimum_portion within the buffered_poly_latlon,
    #then buffer the total_boundaries to the union of those and the original poly
    print('getting sub-jurisdictions for', name)
    admin_lvls = [str(x) for x in range(4,11)] #admin_area below 4 is state/province
    try:
        jurisdictions_latlon = ox.geometries_from_polygon(buffered_poly_latlon, tags={'admin_level':admin_lvls})
    except:
        jurisdictions_latlon = gpd.GeoDataFrame()
    if 'admin_level' in jurisdictions_latlon.columns:
        print(f'found {len(jurisdictions_latlon)} on first pass')
        jurisdictions_utm = jurisdictions_latlon.to_crs(buffered_poly_utm_gdf.crs)
        jurisdictions_utm = gpd.clip(jurisdictions_utm, nearby_land_gdf_utm.unary_union)
        jurisdictions_clipped_utm = jurisdictions_utm.intersection(buffered_poly_utm)
        selection = (jurisdictions_clipped_utm.area / jurisdictions_utm.area) > minimum_portion
        select_jurisdictions_utm = jurisdictions_utm[selection]
        print(f'found {len(select_jurisdictions_utm)} with {minimum_portion} inside area')
        select_jurisdictions_latlon = select_jurisdictions_utm.to_crs(4326)
    else:
        select_jurisdictions_latlon = []
    if len(select_jurisdictions_latlon) > 0:
        total_boundaries_latlon = unary_union([select_jurisdictions_latlon.unary_union, buffered_poly_latlon])
        total_boundaries_utm = unary_union([select_jurisdictions_utm.unary_union, buffered_poly_utm])
    else:
        total_boundaries_latlon = buffered_poly_latlon
        total_boundaries_utm = buffered_poly_utm
        
    #total_boundaries is the area within all admin_areas identified in STEP 1
    #STEP 2: get all jurisdictions within total_boundaries
    #this is so that if we've decided to include a particular admin_area, 
    # we'll also include all admin_areas within it 
    try:
        jurisdictions_latlon = ox.geometries_from_polygon(total_boundaries_latlon, tags={'admin_level':admin_lvls})
    except:
        jurisdictions_latlon = gpd.GeoDataFrame()
        
    if not 'admin_level' in jurisdictions_latlon.columns:
        final_jurisdictions_latlon = []
    else:
        try:
            jurisdictions_latlon = jurisdictions_latlon.loc[('relation',)]
            print(f'found {len(jurisdictions_latlon)} on second pass')
            jurisdictions_utm = jurisdictions_latlon.to_crs(buffered_poly_utm_gdf.crs)
            jurisdictions_utm = gpd.clip(jurisdictions_utm, nearby_land_gdf_utm.unary_union)
            jurisdictions_clipped_utm = jurisdictions_utm.intersection(total_boundaries_utm)
            selection = (jurisdictions_clipped_utm.area / jurisdictions_utm.area) > 0.95
            select_jurisdictions_utm = jurisdictions_utm[selection]
            print(f'found {len(select_jurisdictions_utm)} with 0.95 inside total area')
            selected_levels = []
            for admin_level in select_jurisdictions_utm.admin_level.unique():
                selection = select_jurisdictions_utm[select_jurisdictions_utm.admin_level==admin_level]
                if selection.area.mean() >= level_min_mean_area*1000000:
                    if selection.unary_union.area >= (level_min_coverage * buffered_poly_utm.area):
                        selected_levels.append(admin_level)
                    else:
                        print(f'admin_level={admin_level} excluded: insufficient coverage')
                else:
                    print(f'admin_level={admin_level} excluded: polys too small: avg {selection.area.mean()/1000000}km2')
            final_jurisdictions_utm = select_jurisdictions_utm[select_jurisdictions_utm.admin_level.isin(selected_levels)]
            final_jurisdictions_latlon = final_jurisdictions_utm.to_crs(4326)
            print(f'found {len(final_jurisdictions_latlon)} in acceptable admin levels {selected_levels}')
        except: 
            final_jurisdictions_latlon = []
    
    # get admin_level names, add to dataframe
    level_names_eng = pd.read_csv('input_data/admin_level_names_eng.csv')
    level_names_eng.index = level_names_eng['ISO country code']
    level_names_local = pd.read_csv('input_data/admin_level_names_local.csv')
    level_names_local.index = level_names_local['ISO country code']
    
    if len(final_jurisdictions_latlon) > 0:
        for osmid in final_jurisdictions_latlon.index:
            this_admin_level = final_jurisdictions_latlon.loc[osmid, 'admin_level']
            this_poly = final_jurisdictions_latlon.loc[osmid, 'geometry']
            containers = final_jurisdictions_latlon[
                (final_jurisdictions_latlon.contains(this_poly)) & 
                (final_jurisdictions_latlon.admin_level == this_admin_level)]
            if len(containers) > 1:
                final_jurisdictions_latlon.drop(osmid)
            else:
                try:
                    analysis_areas.loc[new_id,'osmid'] = osmid
                except:
                    import pdb; pdb.set_trace()
                analysis_areas.loc[:,'geometry'].loc[new_id] = final_jurisdictions_latlon.loc[osmid,'geometry']
                #the above hack is necessary because sometimes geometry is a multipolygon
                for attr in ['name','admin_level']:
                    analysis_areas.loc[new_id,attr] = final_jurisdictions_latlon.loc[osmid,attr]
                    analysis_areas.loc[new_id, 'hdc'] = hdc
                level_number = final_jurisdictions_latlon.loc[osmid,'admin_level']
                
                try:
                    level_name_eng = level_names_eng.loc[main_country, f'{level_number}']
                except KeyError:
                    level_name_eng = f'admin_level_{level_number}'
                if type(level_name_eng) != type('string'):
                    level_name_eng = f"admin_level {level_number}"
                try:
                    level_name_local = level_names_local.loc[main_country, f'{level_number}']
                except KeyError:
                    level_name_local = None
                analysis_areas.loc[new_id, 'level_name_eng'] = level_name_eng
                analysis_areas.loc[new_id, 'level_name_local'] = level_name_local
                if type(level_name_local) != type('string'):
                    level_name_full = level_name_eng
                else:
                    level_name_full = f'{level_name_eng} ({level_name_local})'
                analysis_areas.loc[new_id, 'level_name_full'] = level_name_full
                new_id += 1
    
        
    return analysis_areas

def regional_analysis(hdc, 
                      folder_prefix = 'cities_out', 
                      minimum_portion=0.6,
                      prep = True,
                      analyze=True,
                      summarize=True,
                      simplification=0.001, #toposimplification factor
                      cleanup=True, 
                      current_year=2024,
                      ):
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
    analysis_areas = get_jurisdictions(
        hdc, 
        minimum_portion=minimum_portion
        )
    
    analysis_areas.to_file(f'{folder_name}/debug/analysis_areas.gpkg', driver='GPKG')
    
    #3) prepare OSM data files for the agglomeration
    if prep:
        #If we're going to do any access-based indicators
        #Let's make sure to buffer this to include peripheral roads etc for routing
        total_poly_latlon=analysis_areas.unary_union
        total_poly_latlon = shapely.convex_hull(total_poly_latlon)   
        gpd.GeoDataFrame(geometry=[total_poly_latlon], crs=4326).to_file(f'{folder_name}/debug/area_for_osm_extract.gpkg', driver='GPKG')
        prep_from_poly(total_poly_latlon, folder_name, boundary_buffer = 2000)
    
    #4) now actually call the functions and save the results
    if analyze == True:
        geospatial_calctime = pedestriansfirst.spatial_analysis(
            total_poly_latlon,
            hdc,
            analysis_areas.loc[0,'name'],
            folder_name=folder_name,
            )
    else:
        geospatial_calctime = 0
    
    #5) calculate indicator measurement for each analysis area
    if summarize == True:
        start = datetime.datetime.now()
        pedestriansfirst.calculate_indicators(
            analysis_areas, 
            folder_name=folder_name,
            )
        end = datetime.datetime.now()
        analysis_areas.loc[0, 'geospatial_calctime'] = str(geospatial_calctime)
        analysis_areas.loc[0, 'summary_calctime'] = str(end - start)
        
        topo = topojson.Topology(analysis_areas, prequantize=True)
        analysis_areas = topo.toposimplify(0.001).to_gdf()
        analysis_areas.to_file(f'{folder_name}indicator_values.gpkg',driver='GPKG')
        
        nongeospatial_results = analysis_areas.drop('geometry', axis=1, inplace=False)
        nongeospatial_results.to_csv(f'{folder_name}indicator_values.csv')
    
    #clean up big files
    if cleanup:
        for cleanup_filename in ['city.o5m', 'city.pbf','cityhighways.o5m','citywalk.o5m','access/city_ltstagged.pbf']:
            if os.path.exists(f'{folder_name}/temp/{cleanup_filename}'):
                os.remove(f'{folder_name}/temp/{cleanup_filename}')
        for file in os.listdir('cache/'):
            os.remove(f'cache/{file}')
    #import pdb; pdb.set_trace()
        
    
def get_number_jurisdictions(input_folder_prefix = 'cities_out/'):
    """return the total number of jursidictions in all calculated agglomerations
    """
    folders = os.listdir(input_folder_prefix)
    total = 0
    for folder in tqdm(folders):
        data = pd.read_csv(f'{input_folder_prefix}/{folder}/indicator_values.csv')
        if 'admin_level' in data.columns:
            total += len(data[~np.isnan(data.admin_level)])
    print(total)
    return total

def recalculate_blocks(input_folder_prefix = 'cities_out/',current_year=2024,block_patch_length=1000):
    """re-calculate block density 'patches' from existing block geodata
    """
    
    folders = os.listdir(input_folder_prefix)
    for folder in tqdm(folders):
        print(folder)
        blocks_latlon = gpd.read_file(f"{input_folder_prefix}/{folder}/geodata/blocks/blocks_latlon_{current_year}.geojson")
        blocks_utm = ox.project_gdf(blocks_latlon)
        
        block_patches_latlon, block_unbuffered_patches_latlon = pedestriansfirst.make_patches(
            gpd.GeoDataFrame(geometry=[shapely.geometry.box(*blocks_latlon.total_bounds)], crs=4326), 
            blocks_utm.crs, 
            patch_length=block_patch_length
            )
        centroids = blocks_latlon.centroid
        block_unbuf_patches_utm = block_unbuffered_patches_latlon.to_crs(blocks_utm.crs)
        patch_densities = block_unbuffered_patches_latlon
        for patch_idx  in list(patch_densities.index):
            try:
                patch_densities.loc[patch_idx,'block_count'] = centroids.intersects(patch_densities.loc[patch_idx,'geometry']).value_counts()[True]
                #import pdb; pdb.set_trace()
            except KeyError:
                patch_densities.loc[patch_idx,'block_count'] = 0 
        patch_densities_utm = patch_densities.to_crs(blocks_utm.crs)
        patch_densities_utm['density'] = patch_densities_utm.block_count / (patch_densities_utm.area / 1000000)
        patch_densities_latlon = patch_densities_utm.to_crs(epsg=4326)
        try:
            os.rename(f"{input_folder_prefix}/{folder}/geodata/blocks/block_densities_latlon_{current_year}.geojson",
                  f"{input_folder_prefix}/{folder}/geodata/blocks/OLD_block_densities_latlon_{current_year}.geojson")
        except FileNotFoundError:
            pass
        patch_densities_latlon.to_file(f"{input_folder_prefix}/{folder}/geodata/blocks/block_densities_latlon_{current_year}.geojson", driver='GeoJSON')

def make_block_only_folder(input_folder_prefix = 'cities_out/', output_folder_prefix='cities_out_block_data_only/'):
    if not os.path.exists(output_folder_prefix):
        os.makedirs(output_folder_prefix)
    
    folders = os.listdir(input_folder_prefix)
    for folder in tqdm(folders):
        os.makedirs(f'{output_folder_prefix}/{folder}/geodata/blocks')
        shutil.copy(f"{input_folder_prefix}/{folder}/geodata/blocks/block_densities_latlon_2024.geojson",f"{output_folder_prefix}/{folder}/geodata/blocks/block_densities_latlon_2024.geojson")
        shutil.copy(f"{input_folder_prefix}/{folder}/geodata/blocks/blocks_latlon_2024.geojson",f"{output_folder_prefix}/{folder}/geodata/blocks/blocks_latlon_2024.geojson")


def calculate_country_indicators(current_year=2024,
                                 rt_and_pop_years = [1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024, 2025,],
                                 input_folder_prefix = 'cities_out/',
                                 output_folder_prefix = 'countries_out/',
                                 #TODO add years for other indicators with more than one
                                 ):
    """calculate country-level indicators by summarizing city level indicators
    
    This is one of the ugliest functions in the codebase :(
    It can probably be largely reformatted using new (Nov 2024) UCDB data,
    in which agglomerations do NOT cross country borders, and therefore
    much of the complexity won't be needed.
        
    There are different categories of indicator. 
    Some indicators are aggregated to the country level by summing 
    (eg., # of transit stations)
    Others by population-weighted average (eg., People Near or Pop Density)
    
    Also, some are only available for cities with GTFS, or with Rapid Transport
    and other indicators are available for all cities
    
    Note that this includes indicators that are not publicly shown in the 
    Atlas
    """
    
    
    country_bounds = gpd.read_file('input_data/CGAZ/geoBoundaries_ITDPv5_SIMPLIFIED.gpkg')
    
    country_bounds.geometry = country_bounds.make_valid()
    for idx in list(country_bounds.index):
        geometry = country_bounds.loc[idx, 'geometry']
        if geometry.type == 'Polygon':
            country_bounds.loc[idx, 'geometry'] = shapely.geometry.MultiPolygon([geometry])
        elif geometry.type == 'MultiPolygon':
            country_bounds.loc[idx, 'geometry'] = geometry
        elif geometry.type == 'GeometryCollection':
            assert len(list(geometry.geoms)) == 2
            if list(geometry.geoms)[0].type == 'MultiPolygon':
                country_bounds.loc[idx, 'geometry'] = list(geometry.geoms)[0]
            elif list(geometry.geoms)[0].type == 'Polygon':
                country_bounds.loc[idx, 'geometry'] = shapely.geometry.MultiPolygon([list(geometry.geoms)[0]])
            else:
                country_bounds.loc[idx, 'geometry'] = None
                country_bounds.drop(idx, inplace=True)
                print("NNOOOOOOOO")
        else:
            country_bounds.loc[idx, 'geometry'] =None 
            country_bounds.drop(idx, inplace=True)
    
    countries_ISO = list(country_bounds.shapeGroup.unique())
    
    country_regions_organizations = pd.read_csv('input_data/Countries_Regions_Organizations.csv')
    country_regions_organizations.index = country_regions_organizations['ISO Code']
    all_orgs = list(country_regions_organizations.columns)[3:]
    all_regions = list(country_regions_organizations.Region.unique())
    
    if not input_folder_prefix[-1:] == '/':
        input_folder_prefix = input_folder_prefix+'/'
    if not output_folder_prefix[-1:] == '/':
        output_folder_prefix = output_folder_prefix+'/'
    
    #list indicators    
    
    current_year_indicators_avg = [
        'healthcare',
        'schools',
        'h+s',
        'bikeshare',
        'pnab',
        'pnpb',
        'carfree',
        'people_not_near_highways',
        'block_density',
        'pnst',
        ]

    current_year_indicators_sum = [
        'n_points_healthcare',
        'n_points_schools',
        'n_points_bikeshare',
        'highway_km',
        'all_bikeways_km',
        'protected_bikeways_km',
        ]
    
    gtfs_dependent_indicators_avg = [
        'pnft',
        'journey_gap',
        ]
    
    gtfs_dependent_indicators_sum = [
        'n_points_pnft',
        ]
    
    rt_and_pop_indicators_avg = [
        'density',
        'PNrT_all',
        'rtr_all',
        'PNrT_mrt',
        'rtr_mrt',
        'PNrT_lrt',
        'rtr_lrt',
        'PNrT_brt',
        'rtr_brt',
        ]
    
    rt_and_pop_indicators_sum = [
        'total_pop',
        'total_pop_gtfs_cities_only',
        'km_all',
        'stns_all',
        'km_mrt',
        'stns_mrt',
        'km_lrt',
        'stns_lrt',
        'km_brt',
        'stns_brt',
        ]
    
    all_gtfs_dependent = gtfs_dependent_indicators_avg + gtfs_dependent_indicators_sum
    all_currentyear = current_year_indicators_avg + current_year_indicators_sum + gtfs_dependent_indicators_avg + gtfs_dependent_indicators_sum
    all_multiyear = rt_and_pop_indicators_avg + rt_and_pop_indicators_sum
    all_avg = current_year_indicators_avg + gtfs_dependent_indicators_avg + rt_and_pop_indicators_avg
    all_sum = current_year_indicators_sum + gtfs_dependent_indicators_sum + rt_and_pop_indicators_sum
     
    
    full_indicator_names = []
    for indicator in all_currentyear:
        full_indicator_names.append(f'{indicator}_{current_year}')
    for year in rt_and_pop_years:
        for indicator in all_multiyear:
            full_indicator_names.append(f'{indicator}_{year}')
        
            
    #set up dataframes for results
    country_totals = pd.DataFrame(index=countries_ISO, columns=full_indicator_names)
    for country_ISO in countries_ISO:
        name = country_bounds[country_bounds.shapeGroup==country_ISO].shapeName.iloc[0]
        country_totals.loc[country_ISO,'name'] = name
    
    country_totals = country_totals.replace(np.nan,0)
    country_final_values = country_totals.copy()
    region_totals = pd.DataFrame(index = ['world', *all_regions, *all_orgs], columns=full_indicator_names)
    region_totals = region_totals.replace(np.nan,0)
    region_final_values = region_totals.copy()
    
    
    all_cities = gpd.GeoDataFrame(columns=full_indicator_names)
    
    #get data from city-level output
    print('iterating through cities_out/')
    for city_folder in tqdm(os.listdir(f'{input_folder_prefix}')):
        if os.path.exists(f'{input_folder_prefix}{city_folder}/indicator_values.csv'):
            city_results = gpd.read_file(f'{input_folder_prefix}{city_folder}/indicator_values.gpkg')
            #first add to list of full cities
            hdc = city_folder.split('_')[-1]
            all_cities.loc[hdc, 'ID_HDC_G0'] = hdc
            all_cities.loc[hdc, 'name'] = city_results.loc[0,'name']
            all_cities.loc[hdc, 'geometry'] = city_results.loc[0,'geometry']
            for indicator in full_indicator_names:
                if indicator in city_results.columns:
                    all_cities.loc[hdc, indicator] = city_results.loc[0, indicator]
            
            #then calculate by country
            for country in city_results.country.unique():
                if type(country) == type('this is a string, which means it is not np.nan') and country in countries_ISO:
                    if country in country_regions_organizations.index:
                        region = country_regions_organizations.loc[country, 'Region']
                        organizations = list(country_regions_organizations.loc[country,all_orgs][country_regions_organizations.loc[country,all_orgs].notnull()].values)
                        aggregations = ['world',region, *organizations]
                    else:
                        aggregations = ['world']

                    #first total population
                    for year in rt_and_pop_years:
                        total_pop_year = city_results[city_results.country == country][f'total_pop_{year}'].sum()
                        country_totals.loc[country, f'total_pop_{year}'] += total_pop_year
                        for aggregation in aggregations:
                            region_totals.loc[aggregation, f'total_pop_{year}'] += total_pop_year
                        if city_results[city_results.country == country]['has_gtfs'].iloc[0] == True:
                            country_totals.loc[country, f'total_pop_gtfs_cities_only_{year}'] += total_pop_year
                            for aggregation in aggregations:
                                region_totals.loc[aggregation, f'total_pop_gtfs_cities_only_{year}'] += total_pop_year
                    #then indicators based on sums
                    for indicator in full_indicator_names:
                        if indicator[:-5] in all_sum:
                            if not indicator[:9] == 'total_pop':
                                if indicator in city_results.columns:
                                    indicator_total = city_results[city_results.country == country][f'{indicator}'].sum()
                                else:
                                    indicator_total = 0
                                if indicator_total == 'NA':
                                    indicator_total = 0
                                country_totals.loc[country, f'{indicator}'] += indicator_total
                                for aggregation in aggregations:
                                    region_totals.loc[aggregation,f'{indicator}'] += indicator_total
                    #then indicators based on averages
                    for indicator in full_indicator_names:
                        year = indicator[-4:]
                        if indicator in city_results.columns:
                            if indicator[:-5] in all_avg:
                                total_pop_year = city_results[city_results.country == country][f'total_pop_{year}'].sum()                             
                                try:
                                    value = city_results[city_results.country == country][indicator]
                                    multiplied_value = value.astype(float).sum() * total_pop_year
                                    country_totals.loc[country, indicator] += multiplied_value    
                                    for aggregation in aggregations:
                                        region_totals.loc[aggregation,indicator] += multiplied_value 
                                except ValueError:
                                    pass #NA value
        
    #get weighted averages
    print('iterating through countries')
    for country in tqdm(countries_ISO):
        for indicator in full_indicator_names:
            year = indicator[-4:]
            if indicator in country_totals.columns:
                
                if indicator[:-5] in all_sum: #don't need to weight
                    country_final_values.loc[country, indicator] = country_totals.loc[country, indicator]
                
                if indicator[:-5] in all_avg:
                    if indicator[:-5] in gtfs_dependent_indicators_avg:
                        if country_totals.loc[country, f'total_pop_gtfs_cities_only_{year}'] > 0:
                            weighted_avg = country_totals.loc[country, indicator] / country_totals.loc[country, f'total_pop_gtfs_cities_only_{year}']
                        else: 
                            weighted_avg = "n/a"
                    else: #not gtfs-dependent
                        weighted_avg = country_totals.loc[country, indicator] / country_totals.loc[country, f'total_pop_{year}']
                    #import pdb; pdb.set_trace()
                    country_final_values.loc[country, indicator] = weighted_avg
    print('iterating through regions/orgs')
    for region in tqdm(list(region_totals.index)):
        for indicator in full_indicator_names:
            year = indicator[-4:]
            if indicator in region_totals.columns:
                
                if indicator[:-5] in all_sum: #don't need to weight
                    region_final_values.loc[region, indicator] = region_totals.loc[region, indicator]
                
                if indicator[:-5] in all_avg:
                    if indicator[:-5] in gtfs_dependent_indicators_avg:
                        if region_totals.loc[region, f'total_pop_gtfs_cities_only_{year}'] > 0:
                            weighted_avg = region_totals.loc[region, indicator] / region_totals.loc[region, f'total_pop_gtfs_cities_only_{year}']
                        else: 
                            weighted_avg = "n/a"
                    else: #not gtfs-dependent
                        weighted_avg = region_totals.loc[region, indicator] / region_totals.loc[region, f'total_pop_{year}']
                    #import pdb; pdb.set_trace()
                    region_final_values.loc[region, indicator] = weighted_avg
    #save output
    if not os.path.exists(f'{output_folder_prefix}'):
        os.mkdir(f'{output_folder_prefix}')
    country_final_values.to_csv(f'{output_folder_prefix}country_results.csv')
    region_final_values.to_csv(f'{output_folder_prefix}region_results.csv')
    country_geometries = []
    for country in country_final_values.index:
        country_geometries.append(country_bounds[country_bounds.shapeGroup == country].unary_union)
    country_gdf = gpd.GeoDataFrame(country_final_values, geometry=country_geometries, crs=4326)
    country_gdf.to_file(f'{output_folder_prefix}country_results.geojson', driver='GeoJSON')
    all_cities.crs=4326
    all_cities.to_file(f'{output_folder_prefix}all_cities.geojson',driver='GeoJSON')
    pd.DataFrame(all_cities.drop(columns='geometry')).to_csv(f'{output_folder_prefix}all_cities.csv')
    
    #now calculate regional values
    
    
if __name__ == '__main__':
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    ox.utils.config(log_console = False)
    ucdb = gpd.read_file('input_data/ghsl/SMOD_V1s6_opr_P2023_v1_2020_labelUC_DB_release.gpkg')
    ucdb.index =  ucdb['ID_UC_G0']
    hdcs_to_test = [11480 , 461 , 576 , 8265 , 4494] #ITDP "Cycling Cities" below 500k
    hdcs_to_test += list(ucdb[(int(sys.argv[2]) < ucdb['POP_2020'])&(ucdb['POP_2020'] < int(sys.argv[1]))].sort_values('POP_2020', ascending=False).ID_UC_G0)
    for hdc in hdcs_to_test:
        hdc = int(hdc)
        #if len(sys.argv) == 1:
        #divide_by = 1
        #remainder = 0
        # else: 
        divide_by = int(sys.argv[3])
        remainder = int(sys.argv[4])
        print (f"{hdc}%{divide_by}={hdc % divide_by}, compare to {remainder}, {ucdb.loc[hdc,'NAME_MAIN']}, pop {ucdb.loc[hdc,'POP_2020']}")
        if hdc % divide_by == remainder and ucdb.loc[hdc,'NAME_MAIN'] != 'N/A':
            if not os.path.exists(f'cities_out/ghsl_region_{hdc}/indicator_values.csv'):
                if os.path.exists(f'cities_out/ghsl_region_{hdc}/geodata/blocks/blocks_latlon_2024.geojson'):
                    regional_analysis(hdc)#, analyze=False)
                else:
                    regional_analysis(hdc)
    calculate_country_indicators()


