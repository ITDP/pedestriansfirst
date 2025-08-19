#TODO: make this cut out water!
#TODO -- consider customizing this for the USA? an extra buffer or something?
# hdc = 17

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

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

# hdc = "08099" # new york
# hdc = "00507" # tijuana
# hdc = "07532" # chicago
# hdc = "08255" # boston
# hdc = "99999" # kohima
# hdc = "08154" # kohima

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
    # if hdc in [8265, 102]: #kohima, zapopan
    if hdc in ['03105']: #zapopan
        buffer = 40000
    
    #load Urban Centre Database, assign full names
    # here we updated to 2024
    ucdb = gpd.read_file('input_data/ghsl/ghsl_2024.gpkg')
    # ucdb = gpd.read_file('input_data/ghsl/SMOD_V1s6_opr_P2023_v1_2020_labelUC_DB_release.gpkg')
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
    poly_mw_gdf = gpd.GeoDataFrame(geometry=[ghsl_boundaries_mw], crs="4326")
    poly_ll_gdf = poly_mw_gdf.to_crs(4326)
    ghsl_boundaries = poly_ll_gdf.unary_union
    poly_utm_gdf = ox.projection.project_gdf(poly_ll_gdf)
    buffered_poly_utm_gdf = poly_utm_gdf.buffer(buffer)
    buffered_poly_latlon_gdf = buffered_poly_utm_gdf.to_crs(4326)
    buffered_poly_utm = buffered_poly_utm_gdf.unary_union
    buffered_poly_latlon = buffered_poly_latlon_gdf.unary_union
    
    
    # buffered_poly_utm_gdf.to_file("teste_kaue/salvador_boundaires.gpkg")
    
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
    # country_bounds = gpd.read_file("input_data/regions/itdp_boundaries_updated.gpkg")
    
    country_bounds.crs=4326
    country_bounds.geometry = country_bounds.make_valid()
    earth_utm = country_bounds.to_crs(crs = poly_utm_gdf.crs)
    #get land within 10km
    area_for_land_ll = buffered_poly_utm_gdf.buffer(100000).to_crs(4326).unary_union
     # buffered_poly_utm_gdf.to_file("teste_kaue/area_for_land_ll.gpkg")
    
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
    # only filter the country with the largest intersection
    country_overlaps = country_overlaps.loc[country_overlaps['shapeGroup'] == main_country]
    
    # we create analysis areas for "the portion of an agglomeration within each 
    # country
    # May no longer be needed with new UCDB data (coming Nov 2024)
    countries = list(country_overlaps.sort_values('shapeGroup', ascending=True).shapeGroup)
    agglomeration_country_name = ' / '.join(countries)
    analysis_areas.loc[0, 'agglomeration_country_name'] = agglomeration_country_name

    # add country-specific analysis areas
    # I will comment because I don't need this anymore!
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
    ox.settings.overpass_settings = '[out:json][timeout:20]{maxsize}'
    # ox.settings.overpass_settings = '[out:json][timeout:{timeout}]{maxsize}'
    try:
        jurisdictions_latlon = ox.features_from_polygon(buffered_poly_latlon, tags={'admin_level':admin_lvls})
        # remove duplicates
        jurisdictions_latlon = jurisdictions_latlon[~jurisdictions_latlon.index.duplicated(keep='first')]
        # jurisdictions_latlon.to_file('teste_kaue/teste_salvador.gpkg')
    except:
        jurisdictions_latlon = gpd.GeoDataFrame()
    if 'admin_level' in jurisdictions_latlon.columns:
        print(f'found {len(jurisdictions_latlon)} on first pass')
        jurisdictions_utm = jurisdictions_latlon.to_crs(buffered_poly_utm_gdf.crs)
        # UNCOMMENTS LATER
        # jurisdictions_utm = gpd.clip(jurisdictions_utm, nearby_land_gdf_utm.unary_union)
        jurisdictions_clipped_utm = jurisdictions_utm.intersection(buffered_poly_utm)
        selection = (jurisdictions_clipped_utm.area / jurisdictions_utm.area) > minimum_portion
        select_jurisdictions_utm = jurisdictions_utm[selection]
        print(f'found {len(select_jurisdictions_utm)} with {minimum_portion} inside area')
        select_jurisdictions_latlon = select_jurisdictions_utm.to_crs(4326)
        # select_jurisdictions_latlon.to_file('teste_kaue/teste_salvador1.gpkg')
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
    #overpass is used here !!!!!!!!!!!!!!  
        jurisdictions_latlon = ox.features_from_polygon(total_boundaries_latlon, tags={'admin_level':admin_lvls})
        # remove duplicates
        jurisdictions_latlon = jurisdictions_latlon[~jurisdictions_latlon.index.duplicated(keep='first')]
        # teste
        # jurisdictions_latlon1 = ox.features_from_polygon(total_boundaries_latlon, tags={'admin_level':admin_lvls})
        # jurisdictions_latlon.to_file('teste_kaue/teste_salvador2.gpkg')
    except:
        jurisdictions_latlon = gpd.GeoDataFrame()
        
    if not 'admin_level' in jurisdictions_latlon.columns:
        final_jurisdictions_latlon = []
    else:
        try:
            jurisdictions_latlon = jurisdictions_latlon.loc[('relation',)]
            
            print(f'found {len(jurisdictions_latlon)} on second pass')
            jurisdictions_utm = jurisdictions_latlon.to_crs(buffered_poly_utm_gdf.crs)
            # jurisdictions_utm.to_file('teste_kaue/teste_salvadorbefore.gpkg')
            # THE PROBLEM IS ON THIS STEP!!!!
            # jurisdictions_utm = gpd.clip(jurisdictions_utm, nearby_land_gdf_utm.unary_union)
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
            # final_jurisdictions_latlon.to_file('teste_kaue/teste_salvador3.gpkg')
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
          # osmid = final_jurisdictions_latlon.index[310]
            this_admin_level = final_jurisdictions_latlon.loc[osmid, 'admin_level']
            this_poly = final_jurisdictions_latlon.loc[osmid, 'geometry']
            # erro aqui!!!!
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
    
    # analysis_areas.to_csv(f'teste_kaue/teste_{hdc}.csv', sep=',')
    # analysis_areas.to_file(f'teste_kaue/teste_salvador_final.gpkg', driver='GPKG')
    return analysis_areas


# hdcs_to_test = [11480 , 461 , 576 , 8265 , 4494] #ITDP "Cycling Cities" below 500k

# apply function
# teste = get_jurisdictions(hdc = 11480) # this one worked
# teste = get_jurisdictions(hdc = 17)
# teste = get_jurisdictions(hdc = "08154")
# teste.to_file(f'teste_kaue/for_analysis_areas.gpkg', driver='GPKG')
