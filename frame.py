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
    
    
    
def download_ghsl(proj='mw', resolution='100'):
    if not os.path.exists('input_data/'):
        os.mkdir('input_data/')
    if not os.path.exists('input_data/ghsl_data_100m_mw/'):
        os.mkdir('input_data/ghsl_data_100m_mw/')
    for year in tqdm(range(1975, 2031, 5)): 
        letter='E'
        if proj == 'mw':
            name = f'GHS_POP_{letter}{year}_GLOBE_R2023A_54009_{resolution}'
        if not os.path.exists(f'input_data/ghsl_data_100m_mw/{name}/{name}_V1_0.tif'):
            zippath, _ = urllib.request.urlretrieve(f'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/{name}/V1-0/{name}_V1_0.zip')
            with zipfile.ZipFile(zippath, "r") as f:
                f.extractall(f'input_data/ghsl_data_100m_mw/{name}')
    
    
# def from_id_hdc(hdc, folder_prefix = 'cities_out', boundary_buffer = 500, kwargs = {}):
#     if folder_prefix:
#         folder_name = folder_prefix+'/ghsl_'+str(hdc)
#     else:
#         folder_name = str(hdc)
#     poly, name_long, name_short = poly_from_ghsl_hdc(hdc)
#     overpass = prep_from_poly(poly, folder_name, boundary_buffer)
#     kwargs['overpass'] = overpass
#     calctime = pedestriansfirst.spatial_analysis(poly, hdc, name, folder_name, **kwargs)
#     results = pedestriansfirst.calculate_indicators(poly, folder_name)
#     results['calctime'] = calctime
#     print(results)
    
# def from_id_osm(osmid, folder_prefix = 'cities_out', boundary_buffer = 500, kwargs = {}):
#     if folder_prefix:
#         folder_name = folder_prefix+'/osm_'+str(osmid)
#     else:
#         folder_name = str(osmid)
#     poly, name = poly_from_osm_cityid(osmid)
#     overpass = prep_from_poly(poly, folder_name, boundary_buffer)
#     kwargs['overpass'] = overpass
#     calctime = pedestriansfirst.spatial_analysis(poly, osmid, name, folder_name, **kwargs)
#     results = pedestriansfirst.calculate_indicators(poly, folder_name)
#     results['calctime'] = calctime
#     print(results)

# def get_pop_ghsl(city):
#     return city['properties']['P15']

#TODO: make this cut out water!
#TODO -- consider customizing this for the USA? an extra buffer or something?
def get_jurisdictions(hdc,
                      minimum_portion = 0.5, #portion of a jurisdiction that has to be within the poly
                      level_min_mean_area = 2,# min size in km2 for the mean area of a unit at an admin_level
                      level_min_coverage = .0000002, #min coverage of an admin_level of the poly_latlon
                      buffer = 4000, #in m
                      ): 
    
    ucdb = gpd.read_file('input_data/ghsl/SMOD_V1s6_opr_P2023_v1_2020_labelUC_DB_release.gpkg')
    ucdb.index =  ucdb['ID_UC_G0']
    ghsl_boundaries_mw = ucdb.loc[hdc,'geometry']
    name_full = ucdb.loc[hdc,'NAME_LIST']
    all_names = name_full.split('; ')
    name = "The " + " / ".join(all_names[:3]) + ' area'
    if len(name) >= 50:
        name = "The " + " / ".join(all_names[:1]) + ' area'
        if name >= 50:
            name = "The " + all_names[0] + ' area'
    name_short = "The " + all_names[0] + ' area'
    
    poly_mw_gdf = gpd.GeoDataFrame(geometry=[ghsl_boundaries_mw], crs="ESRI:54009")
    poly_ll_gdf = poly_mw_gdf.to_crs(4326)
    ghsl_boundaries = poly_ll_gdf.unary_union
    poly_utm_gdf = ox.project_gdf(poly_ll_gdf)
    buffered_poly_utm_gdf = poly_utm_gdf.buffer(buffer)
    buffered_poly_latlon_gdf = buffered_poly_utm_gdf.to_crs(4326)
    buffered_poly_utm = buffered_poly_utm_gdf.unary_union
    buffered_poly_latlon = buffered_poly_latlon_gdf.unary_union
    
    analysis_areas = gpd.GeoDataFrame()
    new_id = 0
    analysis_areas.loc[new_id, 'name'] = name
    analysis_areas.loc[new_id, 'name_short'] = name_short
    analysis_areas.loc[new_id, 'geometry'] = ghsl_boundaries
    analysis_areas.loc[new_id, 'hdc'] = hdc
    analysis_areas.loc[new_id, 'osmid'] = None
    analysis_areas.loc[new_id, 'level_name'] = 'Agglomeration'
    analysis_areas.crs=4326
    new_id += 1
    
    #now figure out what country it's in
    #and while we're at it, set up country-specific analysis areas
    
    #get CGAZ data here, also use it for clipping coastline
    country_bounds = gpd.read_file('input_data/CGAZ/geoBoundariesCGAZ_ADM0.gpkg')
    country_bounds.crs=4326
    earth_utm = country_bounds.to_crs(crs = poly_utm_gdf.crs)
    #get land within 10km
    area_for_land_ll = buffered_poly_utm_gdf.buffer(100000).to_crs(4326).unary_union
    if country_bounds.intersection(area_for_land_ll).area.sum() >= area_for_land_ll.area * 0.95:
        nearby_land_gdf_utm = buffered_poly_utm_gdf.buffer(100000)
        nearby_land_gdf_ll = buffered_poly_utm_gdf.to_crs(4326)
    else:
        nearby_land_gdf_ll = gpd.clip(country_bounds, area_for_land_ll)
        nearby_land_gdf_utm = nearby_land_gdf_ll.to_crs(buffered_poly_utm_gdf.crs)
            
    country_overlaps = country_bounds.overlay(poly_ll_gdf, how='intersection')
    country_overlaps['land_area'] = country_overlaps.to_crs(poly_utm_gdf.crs).area
    main_country = country_overlaps.sort_values('land_area', ascending=False).iloc[0]['shapeGroup']
    
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
    
        
    
    #First, get all the sub-jusisdictions at least minimum_portion within the buffered_poly_latlon,
    #then buffer the total_boundaries to the union of those and the original poly
    print('getting sub-jurisdictions for', name)
    admin_lvls = [str(x) for x in range(4,11)]
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
        
    #now get all jurisdictions within total_boundaries
    try:
        jurisdictions_latlon = ox.geometries_from_polygon(total_boundaries_latlon, tags={'admin_level':admin_lvls})
    except:
        jurisdictions_latlon = gpd.GeoDataFrame()
        
    if not 'admin_level' in jurisdictions_latlon.columns:
        final_jurisdictions_latlon = []
    else:
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
    
    # get admin_level names, add to dataframe
    level_names_eng = pd.read_csv('input_data/admin_level_names_eng.csv')
    level_names_eng.index = level_names_eng['ISO country code']
    level_names_local = pd.read_csv('input_data/admin_level_names_local.csv')
    level_names_local.index = level_names_local['ISO country code']
    
    if len(final_jurisdictions_latlon) > 0:
        for osmid in final_jurisdictions_latlon.index:
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
                      ):
    
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
    
    #this should happen in get_jurisdictions
    analysis_areas = get_jurisdictions(
        hdc, 
        minimum_portion=minimum_portion
        )
    
    analysis_areas.to_file(f'{folder_name}/debug/analysis_areas.gpkg', driver='GPKG')
    
    if prep:
        #Let's make sure to buffer this to include peripheral roads etc for routing
        total_poly_latlon=analysis_areas.unary_union
        total_poly_latlon = shapely.convex_hull(total_poly_latlon)   
        gpd.GeoDataFrame(geometry=[total_poly_latlon], crs=4326).to_file(f'{folder_name}/debug/area_for_osm_extract.gpkg', driver='GPKG')
        prep_from_poly(total_poly_latlon, folder_name, boundary_buffer = 2000)
    
    #now actually call the functions and save the results
    if analyze == True:
        geospatial_calctime = pedestriansfirst.spatial_analysis(
            total_poly_latlon,
            hdc,
            analysis_areas.loc[0,'name'],
            folder_name=folder_name,
            )
    else:
        geospatial_calctime = 0
    
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
        

def calculate_country_indicators(current_year=2022,
                                 rt_and_pop_years = [1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2025],
                                 input_folder_prefix = 'cities_out/',
                                 output_folder_prefix = 'countries_out/',
                                 #TODO add years for other indicators with more than one
                                 ):
    natural_earth = gpd.read_file('input_data/naturalearth_countries/ne_10m_admin_0_countries.shp')
    countries_ISO = list(natural_earth.ISO_A3.unique())
    
    if not input_folder_prefix[-1:] == '/':
        input_folder_prefix = input_folder_prefix+'/'
    if not output_folder_prefix[-1:] == '/':
        output_folder_prefix = output_folder_prefix+'/'
    
    #list indicators
    current_year_indicators = [
        'healthcare',
        'schools',
        'h+s',
        'bikeshare',
        'pnab',
        'pnpb',
        'carfree',
        'highways',
        'n_points_healthcare',
        'n_points_schools',
        'n_points_bikeshare',
        'people_not_near_highways',
        'highway_km',
        'all_bikeways_km',
        'protected_bikeways_km',
        'block_density',
        ]
    
    gtfs_dependent_indicators = [
        'pnft',
        'n_points_pnft',
        'journey_gap',
        ]
    
    current_year_indicators += gtfs_dependent_indicators
    
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
     
    
    full_indicator_names = []
    for indicator in current_year_indicators:
        full_indicator_names.append(f'{indicator}_{current_year}')
    for year in rt_and_pop_years:
        for indicator in rt_and_pop_indicators_avg:
            full_indicator_names.append(f'{indicator}_{year}')
        for indicator in rt_and_pop_indicators_sum:
            full_indicator_names.append(f'{indicator}_{year}')
        
            
    #set up dataframes for results
    country_totals = pd.DataFrame(index=countries_ISO, columns=full_indicator_names)
    country_totals = country_totals.replace(np.nan,0)
    country_weighted_avgs = country_totals.copy()
    
    all_cities = pd.DataFrame(columns=full_indicator_names)
    
    #get data from city-level output
    print('iterating through cities_out/')
    for city_folder in tqdm(os.listdir(f'{input_folder_prefix}')):
        if os.path.exists(f'{input_folder_prefix}{city_folder}/indicator_values.csv'):
            city_results = pd.read_csv(f'{input_folder_prefix}{city_folder}/indicator_values.csv')
            #first add to list of full cities
            hdc = city_folder.split('_')[-1]
            all_cities.loc[hdc, 'ID_HDC_G0'] = hdc
            all_cities.loc[hdc, 'name'] = city_results.loc[0,'name']
            for indicator in full_indicator_names:
                if indicator in city_results.columns:
                    all_cities.loc[hdc, indicator] = city_results.loc[0, indicator]
            
            #then calculate by country
            for country in city_results.country.unique():
                if type(country) == type('this is a string, which means it is not np.nan'):
                    #indicators based on sums first
                    for year in rt_and_pop_years:
                        total_pop_year = city_results[city_results.country == country][f'total_pop_{year}'].sum()
                        country_totals.loc[country, f'total_pop_{year}'] += total_pop_year
                        if city_results[city_results.country == country]['has_gtfs'].iloc[0] == 'True':
                            country_totals.loc[country, f'total_pop_gtfs_cities_only_{year}'] += total_pop_year
                        for indicator in rt_and_pop_indicators_sum:
                            if not indicator[:9] == 'total_pop':
                                if indicator in city_results.columns:
                                    indicator_total = city_results[city_results.country == country][f'{indicator}_{year}'].sum()
                                else:
                                    indicator_total = 0
                                country_totals.loc[country, f'{indicator}_{year}'] += indicator_total
                    #then the other indicators
                    for indicator in full_indicator_names:
                        year = indicator[-4:]
                        if indicator in city_results.columns:
                            if not indicator[:-5] in rt_and_pop_indicators_sum:
                                total_pop_year = city_results[city_results.country == country][f'total_pop_{year}'].sum()
                                try:
                                    value = city_results[city_results.country == country][indicator].astype(float).sum() * total_pop_year
                                except:
                                    pdb.set_trace()
                                country_totals.loc[country, indicator] += value       
        
    #get weighted averages
    print('iterating through countries')
    for country in tqdm(countries_ISO):
        for indicator in full_indicator_names:
            year = indicator[-4:]
            if indicator in country_totals.columns:
                if indicator[:9] == 'total_pop':
                    country_weighted_avgs.loc[country, indicator] = country_totals.loc[country, indicator]
                elif indicator[:-4] in gtfs_dependent_indicators:
                    weighted_avg = country_totals.loc[country, indicator] / country_totals.loc[country, f'total_pop_gtfs_cities_only_{year}']
                    #import pdb; pdb.set_trace()
                    country_weighted_avgs.loc[country, indicator] = weighted_avg
                else:   
                    weighted_avg = country_totals.loc[country, indicator] / country_totals.loc[country, f'total_pop_{year}']
                    #import pdb; pdb.set_trace()
                    country_weighted_avgs.loc[country, indicator] = weighted_avg    
    #save output
    if not os.path.exists(f'{output_folder_prefix}'):
        os.mkdir(f'{output_folder_prefix}')
    country_weighted_avgs.to_csv(f'{output_folder_prefix}country_results.csv')
    country_geometries = []
    for country in country_weighted_avgs.index:
        country_geometries.append(natural_earth[natural_earth.ISO_A3 == country].unary_union)
    country_gdf = gpd.GeoDataFrame(country_weighted_avgs, geometry=country_geometries, crs=4326)
    country_gdf.to_file(f'{output_folder_prefix}country_results.geojson', driver='GeoJSON')
    all_cities.to_csv(f'{output_folder_prefix}all_cities.csv')
    
            

    
if __name__ == '__main__':
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    ox.utils.config(log_console = False)
    ucdb = gpd.read_file('input_data/ghsl/SMOD_V1s6_opr_P2023_v1_2020_labelUC_DB_release.gpkg')
    ucdb.index =  ucdb['ID_UC_G0']
    #for hdc in ucdb[(int(sys.argv[2]) < ucdb.P15)&(ucdb.P15 < int(sys.argv[1]))].sort_values('P15', ascending=False).ID_HDC_G0:
    for hdc in [
                1031, #quebec
                168, #CDMX
                6312, #delhi
                2269, #toulon
                2038, #dunkirk
                2184, #stavanger
                1397, #flz
                1099, #ciudad este
                24, #SD
                ]:
        hdc = int(hdc)
        #if len(sys.argv) == 1:
        divide_by = 1
        remainder = 0
        # else: 
        #     divide_by = int(sys.argv[3])
        #     remainder = int(sys.argv[4])
        print (f"{hdc}%{divide_by}={hdc % divide_by}, compare to {remainder}, {ucdb.loc[hdc,'NAME_MAIN']}")
        if hdc % divide_by == remainder and ucdb.loc[hdc,'NAME_MAIN'] != 'N/A':
            if not os.path.exists(f'cities_out/ghsl_region_{hdc}/indicator_values.csv'):
                if os.path.exists(f'cities_out/ghsl_region_{hdc}/geodata/blocks/blocks_latlon_2022.geojson'):
                    regional_analysis(hdc)#, analyze=False)
                    calculate_country_indicators()
                else:
                    regional_analysis(hdc)
                    calculate_country_indicators()


