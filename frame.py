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
    ucdb = gpd.read_file('input_data/old_ghsl/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg')
    ucdb.index =  ucdb['ID_HDC_G0']
    boundaries = ucdb.loc[hdc,'geometry']
    name = ucdb.loc[hdc,'UC_NM_MN']
    main_country = ucdb.loc[hdc,'CTR_MN_ISO']
    return boundaries, name, main_country

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
    
    
    
def download_ghsl(proj='mw', resolution='1000'):
    if not os.path.exists('input_data/'):
        os.mkdir('input_data/')
    if not os.path.exists('input_data/ghsl/'):
        os.mkdir('input_data/ghsl/')
    for year in tqdm(range(1975, 2031, 5)): 
        if year <= 2020:
            letter='E'
        else:
            letter='P'
        if proj == 'mw':
            name = f'GHS_POP_{letter}{year}_GLOBE_R2022A_54009_{resolution}'
        else:
            raise ValueError
        zippath, _ = urllib.request.urlretrieve(f'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/{name}/V1-0/{name}_V1_0.zip')
        with zipfile.ZipFile(zippath, "r") as f:
            f.extractall(f'input_data/ghsl/{name}')
    
    
def from_id_hdc(hdc, folder_prefix = 'cities_out', boundary_buffer = 500, kwargs = {}):
    if folder_prefix:
        folder_name = folder_prefix+'/ghsl_'+str(hdc)
    else:
        folder_name = str(hdc)
    poly, name, main_country = poly_from_ghsl_hdc(hdc)
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
def get_jurisdictions(hdc,
                      minimum_portion = 0.6, #portion of a jurisdiction that has to be within the poly
                      level_min_mean_area = 2,# min size in km2 for the mean area of a unit at an admin_level
                      level_min_coverage = .0000002, #min coverage of an admin_level of the poly_latlon
                      buffer = 2000, #in m
                      ): 
    
    ghsl_boundaries, name, main_country = poly_from_ghsl_hdc(hdc)
    
    poly_utm_gdf = ox.projection.project_gdf(gpd.GeoDataFrame(geometry=[ghsl_boundaries], crs=4326))
    buffered_poly_utm_gdf = poly_utm_gdf.buffer(buffer)
    buffered_poly_latlon_gdf = buffered_poly_utm_gdf.to_crs(4326)
    buffered_poly_utm = buffered_poly_utm_gdf.unary_union
    buffered_poly_latlon = buffered_poly_latlon_gdf.unary_union
    
    analysis_areas = gpd.GeoDataFrame()
    new_id = 0
    analysis_areas.loc[new_id,'name'] = name
    analysis_areas.loc[new_id, 'geometry'] = ghsl_boundaries
    analysis_areas.loc[new_id, 'hdc'] = hdc
    analysis_areas.loc[new_id, 'osmid'] = None
    analysis_areas.loc[new_id, 'level_name'] = 'Agglomeration'
    analysis_areas.crs=4326
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
            analysis_areas.loc[new_id,'name'] = area[1].name_muni
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
        
    #get natural_earth data here, both to use it for clipping coastline
    #and for countries later
    natural_earth = gpd.read_file('input_data/naturalearth_countries/ne_10m_admin_0_countries.shp')
    earth_utm = natural_earth.to_crs(crs = poly_utm_gdf.crs)
    #get land within 10km
    area_for_land_ll = buffered_poly_utm_gdf.buffer(100000).to_crs(4326).unary_union
    if natural_earth.intersection(area_for_land_ll).area.sum() >= area_for_land_ll.area * 0.95:
        nearby_land_gdf_utm = buffered_poly_utm_gdf.buffer(100000)
        nearby_land_gdf_ll = buffered_poly_utm_gdf.to_crs(4326)
    else:
        nearby_land_gdf_ll = gpd.clip(natural_earth, area_for_land_ll)
        nearby_land_gdf_utm = nearby_land_gdf_ll.to_crs(buffered_poly_utm_gdf.crs)
        
    
    #First, get all the sub-jusisdictions at least minimum_portion within the buffered_poly_latlon,
    #then buffer the total_boundaries to the union of those and the original poly
    print('getting sub-jurisdictions')
    admin_lvls = [str(x) for x in range(4,11)]
    try:
        jurisdictions_latlon = ox.geometries_from_polygon(buffered_poly_latlon, tags={'admin_level':admin_lvls})
    except:
        jurisdictions_latlon = gpd.GeoDataFrame()
    if 'admin_level' in jurisdictions_latlon.columns:
        #jurisdictions_latlon = jurisdictions_latlon.loc[('relation',)]
        #I forget what that was doing
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
    jurisdictions_latlon = ox.geometries_from_polygon(total_boundaries_latlon, tags={'admin_level':admin_lvls})
    if not 'admin_level' in jurisdictions_latlon.columns:
        final_jurisdictions_latlon = []
    else:
        #jurisdictions_latlon = jurisdictions_latlon.loc[('relation',)]
        #I forget what that was doing
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
                analysis_areas.loc[new_id,'osmid'] = osmid[1]
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
                level_name_full = f'{level_name_eng} ({level_name_local}'
            analysis_areas.loc[new_id, 'level_name_full'] = level_name_full
            new_id += 1
            
    country_overlaps = natural_earth.overlay(gpd.GeoDataFrame(geometry=[ghsl_boundaries], crs=4326), how='intersection')
    for idx in country_overlaps.index:
        analysis_areas.loc[new_id,'country'] = country_overlaps.loc[idx,'ISO_A3']
        analysis_areas.loc[:,'geometry'].loc[new_id] = country_overlaps.loc[idx,'geometry']
        new_id += 1
        
    return analysis_areas

def regional_analysis(hdc, 
                      folder_prefix = 'cities_out', 
                      minimum_portion=0.6,
                      analyze=True,
                      summarize=True,
                      simplification=0.001, #toposimplification factor
                      cleanup=False,
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
    #import pdb; pdb.set_trace()
        

#all cities in descending order
#TODO: make use gpd instead of fiona lol
# def all_cities():
#     with fiona.open('GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp','r') as ucdb:
#         cities = list(ucdb)
#     cities.sort(key=get_pop_ghsl, reverse = True)
#     for city in cities:
#         if os.path.exists('all_results.json'):
#             with open('all_results.json','r') as in_file:
#                 all_results = json.load(in_file)
#         else:
#             all_results = {}
#         if not str(city['properties']['ID_HDC_G0']) in all_results.keys():
#             if not str(city['properties']['ID_HDC_G0']) == '4541': #there's one city in south sudan that doesn't work right.
#                 results = from_id_hdc(city['properties']['ID_HDC_G0'])
#                 all_results.update({city['properties']['ID_HDC_G0']:results})
#                 with open('all_results.json','w') as out_file:
#                     json.dump(all_results, out_file)
            
def calculate_country_indicators(current_year=2022,
                                 rt_and_pop_years = [1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2025],
                                 #TODO add years for other indicators with more than one
                                 ):
    natural_earth = gpd.read_file('input_data/naturalearth_countries/ne_10m_admin_0_countries.shp')
    countries_ISO = list(natural_earth.ISO_A3.unique())
    
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
    
    rt_and_pop_indicators = [
        'total_pop',
        'total_pop_gtfs_cities_only',
        'density',
        'PNrT_all',
        'km_all',
        'stns_all',
        'rtr_all',
        'PNrT_mrt',
        'km_mrt'
        'stns_mrt',
        'rtr_mrt',
        'PNrT_lrt',
        'km_lrt'
        'stns_lrt',
        'rtr_lrt',
        'PNrT_brt',
        'km_brt',
        'stns_brt',
        'rtr_brt',
        ]
     
    
    full_indicator_names = []
    for indicator in current_year_indicators:
        full_indicator_names.append(f'{indicator}_{current_year}')
    for year in rt_and_pop_years:
        for indicator in rt_and_pop_indicators:
            full_indicator_names.append(f'{indicator}_{year}')
            
    #set up dataframes for results
    country_totals = pd.DataFrame(index=countries_ISO, columns=full_indicator_names)
    country_totals = country_totals.replace(np.nan,0)
    country_weighted_avgs = country_totals.copy()
    
    #get data from city-level output
    print('iterating through cities_out/')
    for city_folder in tqdm(os.listdir('cities_out/')):
        if os.path.exists(f'cities_out/{city_folder}/indicator_values.csv'):
            city_results = pd.read_csv(f'cities_out/{city_folder}/indicator_values.csv')
            for country in city_results.country.unique():
                if type(country) == type('this is a string, which means it is not np.nan'):
                    #total pop first
                    for year in rt_and_pop_years:
                        total_pop_year = city_results[city_results.country == country][f'total_pop_{year}'].sum()
                        country_totals.loc[country, f'total_pop_{year}'] += total_pop_year
                        if city_results[city_results.country == country]['has_gtfs'].iloc[0] == 'True':
                            country_totals.loc[country, f'total_pop_gtfs_cities_only_{year}'] += total_pop_year
                    #then the other indicators
                    for indicator in full_indicator_names:
                        year = indicator[-4:]
                        if indicator in city_results.columns:
                            if not indicator[:9] == 'total_pop':
                                total_pop_year = country_totals.loc[country, f'total_pop_{year}']
                                value = city_results[city_results.country == country][indicator].sum() * total_pop_year
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
    if not os.path.exists('country_results/'):
        os.mkdir('country_results')
    country_weighted_avgs.to_csv('country_results/country_results.csv')
    country_geometries = []
    for country in country_weighted_avgs.index:
        country_geometries.append(natural_earth[natural_earth.ISO_A3 == country].unary_union)
    country_gdf = gpd.GeoDataFrame(country_weighted_avgs, geometry=country_geometries, crs=4326)
    country_gdf.to_file('country_results/country_results.geojson', driver='GeoJSON')
    
            

    
if __name__ == '__main__':
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    ox.utils.config(log_console = False)
    
    
    
    test_cities = [

        #2185,
        #11640,
        #13043,
        
        
       #1210,
       #1303,
       #1420,
       #21,
       #2851,
       #855,
       #3751,
       #9691,
       #855,
        
        
#     855	,
#                 11498	,
# 57	,
# 88	,
# 561	,
# 4351	,
# 4309	,
# 4359	,
# 12696	,
# 12080	,
# 2051	,
# 3541	,
# 1361	,
# 1445	,
# 1406	,
# 154	,
# 200	,
# 21	,
# 350	,
# 1105	,
# 931	,
# 5134	,
# 4172	,
# 3902	,
# 4427	,
# 4608	,
# 9691	,
# 7041	,
# 10076	,
# 8050	,
# 6522	,
# 11862	,
# 1709	,
# 2806	,
# 2980	,
# 2559	,
# 13039	,
# 3562	,
# 1372	,
# 8675	,
# 7041	,
# 6845	,
# 12394	,
# 12384	,
# 1022	,
# 1009	,
# 828	,
# 2749	,
# 10687	,
# 11223	,
# 14	,
        ]
    
    for cityid in test_cities:
        regional_analysis(cityid)

    ucdb = gpd.read_file('input_data/old_ghsl/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg')
    for hdc in ucdb[(500000 < ucdb.P15)&(ucdb.P15 < 1000000)].sort_values('P15', ascending=False).ID_HDC_G0:
        hdc = int(hdc)
        if not os.path.exists(f'cities_out/ghsl_region_{hdc}/indicator_values.csv'):
            if os.path.exists(f'cities_out/ghsl_region_{hdc}/geodata/blocks/blocks_latlon_2022.geojson'):
                regional_analysis(hdc, analyze=False)
            else:
                regional_analysis(hdc)
        calculate_country_indicators()


