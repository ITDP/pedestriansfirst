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



def calculate_country_indicators(current_year=2024,
                                 rt_and_pop_years = [1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2023, 2024, 2025],
                                 input_folder_prefix = '/media/kauebraga/data/pedestriansfirst/cities_out/',
                                 output_folder_prefix = '/media/kauebraga/data/pedestriansfirst/countries_out/',
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
    # un_countries =      gpd.read_file("input_data/regions/un_countries.gpkg")
    # un_continents =     gpd.read_file("input_data/regions/un_continents.gpkg")
    # un_regions =        gpd.read_file("input_data/regions/un_regions.gpkg")
    # worldbank_regions = gpd.read_file("input_data/regions/worldbank_regions.gpkg")
    # worldbank_income =  gpd.read_file("input_data/regions/worldbank_income.gpkg")
    
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
    
    
    # open the table with the relationshuip with the countries and the other regions
    country_regions_organizations = pd.read_csv("input_data/regions/regions_lookup.csv")
    country_regions_organizations.index = country_regions_organizations['ISO_Alpha3_Code']
    country_regions_organizations["UN_Continent"] = 'UN Continents / ' + country_regions_organizations["UN_Continent"]
    country_regions_organizations["UN_Region"] = 'UN Regions / ' + country_regions_organizations["UN_Region"]
    country_regions_organizations["World_Bank_Region"] = 'World Bank Regions / ' + country_regions_organizations["World_Bank_Region"]
    country_regions_organizations["World_Bank_Income_Group"] = 'World Bank Income / ' + country_regions_organizations["World_Bank_Income_Group"]
    all_un_continents = list(country_regions_organizations.UN_Continent.unique())
    all_un_regions = list(country_regions_organizations.UN_Region.unique())
    all_worldbank_regions = list(country_regions_organizations.World_Bank_Region.unique())
    all_worldbank_income = list(country_regions_organizations.World_Bank_Income_Group.unique())
    # country_regions_organizations = pd.read_csv('input_data/Countries_Regions_Organizations.csv')
    # country_regions_organizations.index = country_regions_organizations['ISO Code']
    # all_orgs = list(country_regions_organizations.columns)[3:]
    # all_regions = list(country_regions_organizations.Region.unique())
    
    if not input_folder_prefix[-1:] == '/':
        input_folder_prefix = input_folder_prefix+'/'
    if not output_folder_prefix[-1:] == '/':
        output_folder_prefix = output_folder_prefix+'/'
    
    #list indicators    
    
    current_year_indicators_avg = [
        'healthcare',
        'schools',
        'hs',
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
    
    # region_totals = pd.DataFrame(index = ['world', *all_regions, *all_orgs], columns=full_indicator_names)
    region_totals = pd.DataFrame(index = ['world', *all_un_continents, *all_un_regions, *all_worldbank_regions, *all_worldbank_income], columns=full_indicator_names)
    region_totals = region_totals.replace(np.nan,0)
    region_final_values = region_totals.copy()
    
    
    all_cities = gpd.GeoDataFrame(columns=full_indicator_names)
    
    #get data from city-level output
    print('iterating through cities_out/')
    for city_folder in tqdm(os.listdir(f'{input_folder_prefix}')):
        if os.path.exists(f'{input_folder_prefix}{city_folder}/indicator_values_{current_year}.csv'):
            # city_folder = "ghsl_region_02732"
            city_results = gpd.read_file(f'{input_folder_prefix}{city_folder}/indicator_values_{current_year}.gpkg')
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
                # country = "JPN"
                if type(country) == type('this is a string, which means it is not np.nan') and country in countries_ISO:
                    if country in country_regions_organizations.index:
                        # region = country_regions_organizations.loc[country, 'Region']
                        # organizations = list(country_regions_organizations.loc[country,all_orgs][country_regions_organizations.loc[country,all_orgs].notnull()].values)
                        un_continent = country_regions_organizations.loc[country, 'UN_Continent']
                        un_region = country_regions_organizations.loc[country, 'UN_Region']
                        worldbank_region = country_regions_organizations.loc[country, 'World_Bank_Region']
                        worldbank_income = country_regions_organizations.loc[country, 'World_Bank_Income_Group']
                        aggregations = ['world',un_continent, un_region, worldbank_region, worldbank_income]
                    else:
                        aggregations = ['world']

                    #first total population
                    for year in rt_and_pop_years:
                        # year = 2025
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
        # country = "BRA"
        for indicator in full_indicator_names:
            # indicator = "schools_2023"
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
        # region = "Northern America"
        for indicator in full_indicator_names:
            # indicator = "pnft_2023"
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
    country_final_values.to_csv(f'{output_folder_prefix}country_results_{current_year}.csv')
    region_final_values.to_csv(f'{output_folder_prefix}region_results_{current_year}.csv')
    country_geometries = []
    for country in country_final_values.index:
        country_geometries.append(country_bounds[country_bounds.shapeGroup == country].unary_union)
    country_gdf = gpd.GeoDataFrame(country_final_values, geometry=country_geometries, crs=4326)
    country_gdf.to_file(f'{output_folder_prefix}country_results_{current_year}.geojson', driver='GeoJSON')
    all_cities.crs=4326
    all_cities.to_file(f'{output_folder_prefix}all_cities_{current_year}.geojson',driver='GeoJSON')
    pd.DataFrame(all_cities.drop(columns='geometry')).to_csv(f'{output_folder_prefix}all_cities_{current_year}.csv')
