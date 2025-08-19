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
# hdc = 8154 # fortaleza
# hdc = 8099 # ny

def calculate_isochrones(hdc,
                          distances = { #network buffers, in meters
                            'healthcare': 1000,
                            'schools': 1000,
                            'libraries': 1000,
                            'bikeshare': 300,
                            'pnft': 500,
                            'pnrt': 1000,
                            'special': 250,
                            'pnpb': 250,
                            'pnab': 250,
                            'highways':500,
                            },
                          to_test = [
                           'healthcare',
                           'schools',
                           'h+s',
                           #'libraries',
                           'bikeshare',
                           # 'carfree',
                           'blocks',
                           'density',
                           # 'pnft', # VOLTAR PRO PADRAO!!!
                           # 'pnrt',
                           'pnpb', #protected bikeways
                           'pnab', #all bikeways
                           # 'pnst', #combo transit + bike
                           'highways',
                           #'journey_gap',
                           #'access',
                           #'transport_performance',
                           #'special',
                           ],
                            buffer_dist=100,#buffer out from nodes, m)        
                            boundary_buffer = 1000, #m
                            services_simplification = 10, #m TODO replace with gpd simplification?
                                                  current_year = 2024,
                      patch_length = 40000, #m
                      block_patch_length = 1000, #m
                      boundary_buffer = 1000, #m
                      blocks_simplification = 0.0001, #topo-simplification
                      services_simplification = 10, #m TODO replace with gpd simplification?
                      access_resolution = 2000, #m
                      
  folder_name = 'cities_out'+'/ghsl_region_'+str(hdc)+'/'
                      
  useful_tags = ox.settings.useful_tags_way + ['cycleway', 'cycleway:left', 'cycleway:right', 'cycleway:both', 'bicycle']
  ox.config(use_cache=True, log_console=True, useful_tags_way=useful_tags)
      
    # open boundaries
    boundaries = gpd.read_file(f'{folder_name}/debug/area_for_osm_extract.gpkg')
    
    boundaries_latlon = gpd.GeoDataFrame(geometry = boundaries['geometry'].values)
    boundaries_latlon.crs = {'init':'epsg:4326'}
    longitude = round(numpy.mean(boundaries_latlon.geometry.centroid.x),10)
    utm_zone = int(math.floor((longitude + 180) / 6) + 1)
    utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
    boundaries_utm = boundaries_latlon.to_crs(utm_crs)
    boundaries_utm.geometry = boundaries_utm.geometry.buffer(boundary_buffer)
    boundaries_latlon = boundaries_utm.to_crs(epsg=4326)
    boundaries_mw = boundaries_utm.to_crs("ESRI:54009")
    boundaries = boundaries_latlon.geometry.unary_union
    
    if 'h+s' in to_test:
    for part in ['healthcare', 'schools']:
        if not part in to_test:
            to_test.append(part)
                            
  
  # open the services
  service_point_locations = {'schools' : gpd.read_file(f'frame1/new_funs/data/schools.gpkg'),
                            'healthcare' : gpd.read_file(f'frame1/new_funs/data/healthcare.gpkg'),
                            'libraries' : gpd.read_file(f'frame1/new_funs/data/libraries.gpkg'),
                            'schools' : gpd.read_file(f'frame1/new_funs/data/schools.gpkg')}
  if pnft 
    service_point_locations = {'schools' : gpd.read_file(f'frame1/new_funs/data/schools.gpkg'),
                              'healthcare' : gpd.read_file(f'frame1/new_funs/data/healthcare.gpkg'),
                              'libraries' : gpd.read_file(f'frame1/new_funs/data/libraries.gpkg'),
                              'schools' : gpd.read_file(f'frame1/new_funs/data/schools.gpkg'),
                              'pnft' : gpd.read_file(f'frame1/new_funs/data/pnft.gpkg')}
    
                            
                
  if 'highways' in to_test:
    all_highway_lines = []
    all_highway_areas = []
    
    
        
  if 'blocks' in to_test:
      outblocks = []
      block_counts = []
  
  quilt_isochrone_polys = {}
  for service in to_test:
      quilt_isochrone_polys[service] = False
      if 'pnab' in to_test or 'pnpb' in to_test:
          quilt_isochrone_polys['pnab'] = False
          quilt_isochrone_polys['pnpb'] = False
          quilt_allbike = gpd.GeoDataFrame(geometry=[], crs=utm_crs)
          quilt_protectedbike = gpd.GeoDataFrame(geometry=[], crs=utm_crs)
    
    
  # create the network
      
    path = folder_name + 'temp/cityhighways.osm'
    G_allroads = ox.graph_from_xml(path, simplify=False, retain_all=True)
    #label links that are not in tunnels
    #so highway identification works properly later
    #
    for x in G_allroads.edges:
        if 'tunnel' not in G_allroads.edges[x].keys():
            G_allroads.edges[x]['tunnel'] = "no" 
            
    #then simplify
    G_allroads = ox.simplify_graph(G_allroads)
    G_allroads.remove_nodes_from(list(nx.isolates(G_allroads)))

    
    # project                
    if len(G_allroads.edges) > 0 and len(G_allroads.nodes) > 0:
        G_allroads = ox.project_graph(G_allroads, to_crs=utm_crs)
        
    # transform the network to gdfs: will be neccesary for process_bike and get_highways
    nodes, edges = ox.graph_to_gdfs(G_allroads)
        
  # process pnab/pnpb
   center_nodes = {}
    if 'pnab' in to_test or 'pnpb' in to_test:
      process_bike(edges, quilt_allbike, quilt_protectedbike, center_nodes) 
  
  # process highways
    if 'highways' in to_test:
    # SLOW HERE !!!!!!!!!!!!
      highway_lines_gdf_ll = get_highways(nodes, edges)
      # highway_lines_gdf_ll = get_service_locations.get_highways(G_allroads)
      if highway_lines_gdf_ll is not None and len(list(highway_lines_gdf_ll.geometry)) > 0:
          all_highway_lines += list(highway_lines_gdf_ll.geometry)
          highway_lines_utm = ox.project_gdf(highway_lines_gdf_ll)
          highway_areas_utm = highway_lines_utm.buffer(distances['highways'])
          highway_areas_ll = highway_areas_utm.to_crs(4326)
          all_highway_areas +=list(highway_areas_ll.geometry)
  
  # identify the nodes from each service
  
  for service in service_point_locations.keys():
        if service in ['healthcare','schools','libraries','bikeshare','pnft','special']:
          # crop the points (schools hospitals etc) to the patch boundaries?
            patch_points = service_point_locations[service]
            patch_points.crs=4326
            patch_points_utm = patch_points.to_crs(utm_crs)
            patch_points_utm = patch_points_utm[(~patch_points_utm.is_empty)
                                                &(patch_points_utm.geometry!=None)]
            
            center_nodes[service]  = ox.distance.nearest_nodes(
                G_allroads, 
                patch_points_utm.geometry.x, 
                patch_points_utm.geometry.y,
                return_dist=False)
                
                
                
  # Get polygons - isocrhones?
  
      
  testing_services = []
  for service in ['healthcare', 'schools', 'libraries', 'bikeshare']:
      if service in to_test:
          testing_services.append(service)
          
  isochrone_polys = {}
  testing_services.append('pnpb')
  testing_services.append('pnab')
  for service in testing_services:
      if service in center_nodes.keys():
          # service = 'pnab'
          print(f'getting polygons for {service}, {len(center_nodes[service])} center_nodes')
          isochrone_polys[service] = proper_iso_polys(
          # isochrone_polys[service] = isochrones.proper_iso_polys(
              G_allroads, 
              center_nodes[service],
              distance=distances[service],                                           
              buffer=buffer_dist, 
              infill=2500)    
              
  for service in isochrone_polys.keys():
    if service not in quilt_isochrone_polys.keys() or not quilt_isochrone_polys[service]:
        quilt_isochrone_polys[service] = isochrone_polys[service]
    elif isochrone_polys[service]:
        quilt_isochrone_polys[service] = shapely.ops.unary_union([quilt_isochrone_polys[service],isochrone_polys[service]])

  # get polygons for rapid transit
  if 'pnrt' in to_test:
      for stn_idx in rt_stns.index:
          if unbuffered_patch.contains(rt_stns.loc[stn_idx, 'geometry']):
              stn_utm = rt_stns_utm.loc[stn_idx, 'geometry']
              center_node = ox.distance.nearest_nodes(
                  G_allroads, 
                  stn_utm.x, 
                  stn_utm.y,
                  return_dist=False)
              iso_poly = proper_iso_polys(
                  G_allroads, 
                  [center_node],
                  distance=distances['pnrt'],                                           
                  buffer=buffer_dist, 
                  infill=2500)
              rt_isochrones_utm.loc[stn_idx,'geometry'] = iso_poly
  
  if 'blocks' in to_test:
      if G_allroads and len(G_allroads.edges) > 0:
          
          # streets = ox.utils_graph.graph_to_gdfs(G_allroads, nodes = False)
          streets = edges
          streets = shapely.geometry.MultiLineString(list(streets.geometry))
          merged = shapely.ops.linemerge(streets)
          if merged:
              borders = shapely.ops.unary_union(merged)
              blocks = list(shapely.ops.polygonize(borders))
              all_blocks = []
              selected_areas = []
              for block in blocks:
                  if 500 < block.area: #< 200000000:
                      if block.interiors:
                          block = shapely.geometry.Polygon(block.exterior)
                      if block.centroid.within(unbuffered_patch_utm):
                          area = round(block.area, 3)
                          perim = round(block.length, 3)
                          lemgth = round((perim * perim) / area, 3)
                          if blocks_simplification:
                              block = block.simplify(blocks_simplification)
                          all_blocks.append((block, area, perim, lemgth))
                          if (lemgth > 75) and (1000 > area or area > 1000000): #obvious outliers
                              pass
                          elif (lemgth > 30 ) and (3000 > area):
                              pass
                          else:
                              selected_areas.append(area)
              outblocks += all_blocks
              print(f'cut {len(all_blocks)} blocks')
              block_counts.append(len(all_blocks))
          else:
              block_counts.append(0)
              print('not merged!')
      else:
          block_counts.append(0)
          
          
          
          
          
          
  del G_allroads
  gc.collect()
  
  debugcounter = 1
  print(debugcounter); debugcounter+=1
  
  
  
  
  
  #start saving files
  geodata_subfolders = []
  for service in testing_services:
      geodata_subfolders.append(service+'/')
      geodata_subfolders.append(service+'_points/')
  geodata_subfolders += [
      'h+s',
      'protectedbike/',
      'allbike/',
      'carfree/',
      'blocks/',
      'allhwys/',
      'buffered_hwys',
      ]
  for subfolder in geodata_subfolders:
      if not os.path.exists(f"{folder_name}geodata/{subfolder}/"):
          os.mkdir(f"{folder_name}geodata/{subfolder}/")
  
  for service in testing_services:
      
      print(service)
      if quilt_isochrone_polys[service]:
          service_utm = gpd.GeoDataFrame(geometry = [quilt_isochrone_polys[service]],
                                         crs=utm_crs)
          #TODO: check validity / all polygons?
          service_utm.geometry = service_utm.geometry.simplify(services_simplification)
          service_utm = gpd.overlay(service_utm ,boundaries_utm, how='intersection')
          service_latlon = service_utm.to_crs(epsg=4326)
          service_latlon.to_file(f"{folder_name}geodata/{service}/{service}_latlon_{current_year}.geojson", driver='GeoJSON')
      if service in service_point_locations.keys():
          service_point_locations[service].to_file(f"{folder_name}geodata/{service}_points/{service}_points_latlon_{current_year}.geojson", driver='GeoJSON')
        
  
  if 'pnpb' in to_test or 'pnab' in to_test:
      if not quilt_protectedbike.empty:
          quilt_protectedbike = quilt_protectedbike.to_crs(4326)
          merged_protectedbike = quilt_protectedbike.intersection(boundaries)
          merged_protectedbike = gpd.GeoDataFrame(geometry = [merged_protectedbike.unary_union], crs=4326)
          merged_protectedbike.to_file(f"{folder_name}geodata/protectedbike/protectedbike_latlon_{current_year}.geojson",driver='GeoJSON')
      if not quilt_allbike.empty:
          quilt_allbike = quilt_allbike.to_crs(4326)
          merged_allbike = quilt_allbike.intersection(boundaries)
          merged_allbike = gpd.GeoDataFrame(geometry = [merged_allbike.unary_union], crs=4326)
          merged_allbike.to_file(f"{folder_name}geodata/allbike/allbike_latlon_{current_year}.geojson",driver='GeoJSON')

      
  if 'highways' in to_test:
      all_hwys_multiline = shapely.ops.unary_union(all_highway_lines)
      all_hwys_latlon = gpd.GeoDataFrame(geometry=[all_hwys_multiline], crs=4326)
      all_hwys_latlon.to_file(f"{folder_name}geodata/allhwys/allhwys_latlon_{current_year}.geojson",driver='GeoJSON')
      
      all_hwys_multipoly = shapely.ops.unary_union(all_highway_areas)
      all_hwy_buffs_latlon = gpd.GeoDataFrame(geometry=[all_hwys_multipoly], crs=4326)
      all_hwy_buffs_latlon.to_file(f"{folder_name}geodata/buffered_hwys/buffered_hwys_latlon_{current_year}.geojson",driver='GeoJSON')
      
      del all_hwys_latlon
      del all_hwys_multiline
      del all_hwys_multipoly
      del all_hwy_buffs_latlon
      
      gc.collect()
      
  if 'carfree' in to_test:
      carfree_latlon = gpd.GeoDataFrame(geometry = citywide_carfree)
      #just a latlon list of points
      carfree_latlon.crs = {'init':'epsg:4326'}
      carfree_utm = carfree_latlon.to_crs(utm_crs)
      carfree_utm.geometry = carfree_utm.geometry.buffer(100)
      #this is the analysis, the 100m buffer
      carfree_utm = gpd.GeoDataFrame(geometry = [shapely.ops.unary_union(carfree_utm.geometry)], crs=utm_crs)
      carfree_utm.geometry = carfree_utm.geometry.simplify(services_simplification)
      carfree_utm = gpd.overlay(carfree_utm ,boundaries_utm, how='intersection')
      carfree_latlon = carfree_utm.to_crs('epsg:4326')
      carfree_latlon.to_file(f"{folder_name}geodata/carfree/carfree_latlon_{current_year}.geojson",driver='GeoJSON')
     
      
  if 'h+s' in to_test:
      if quilt_isochrone_polys['healthcare'] and quilt_isochrone_polys['schools']:
          service = 'h+s'
          intersect = shapely.ops.unary_union([quilt_isochrone_polys['healthcare'].intersection(quilt_isochrone_polys['schools'])])
          if type(intersect) == shapely.geometry.collection.GeometryCollection:
              if not intersect.is_empty:
                  try:
                      intersect = [obj for obj in intersect if type(obj) == shapely.geometry.polygon.Polygon]
                  except TypeError: #intersect is a GeometryCollection
                      intersect = [obj for obj in intersect.geoms if type(obj) == shapely.geometry.polygon.Polygon]
                  intersect = shapely.geometry.MultiPolygon(intersect)
          hs_utm = gpd.GeoDataFrame(geometry = [intersect], crs=utm_crs)
          if hs_utm.geometry.area.sum() != 0:
              hs_utm = gpd.overlay(hs_utm ,boundaries_utm, how='intersection')
              hs_utm.geometry = hs_utm.geometry.simplify(services_simplification)
              hs_latlon = hs_utm.to_crs(epsg=4326)
              hs_latlon.to_file(f"{folder_name}geodata/{service}/{service}_latlon_{current_year}.geojson", driver='GeoJSON')
  
  if 'pnrt' in to_test:
      if not os.path.exists(folder_name+'geodata/rapid_transit/'):
          os.mkdir(folder_name+'geodata/rapid_transit/')
      if 'rt_mode' in rt_isochrones_utm.columns:
          rt_isochrones_latlon = rt_isochrones_utm.to_crs(4326)
          for year in years:
              if not os.path.exists(f'{folder_name}geodata/rapid_transit/{year}/'):
                  os.mkdir(f'{folder_name}geodata/rapid_transit/{year}/')
              #this could probably be more elegant
              mode_selectors = {
                  'brt_atgrade': rt_isochrones_latlon['rt_mode'] == 'brt_atgrade',
                  'brt_gradesep': rt_isochrones_latlon['rt_mode'] == 'brt_gradesep',
                  'brt': rt_isochrones_latlon['rt_mode'].isin(['brt_atgrade','brt_gradesep']),
                  'lrt_atgrade': rt_isochrones_latlon['rt_mode'] == 'lrt_atgrade',
                  'lrt_gradesep': rt_isochrones_latlon['rt_mode'] == 'lrt_gradesep',
                  'lrt': rt_isochrones_latlon['rt_mode'].isin(['lrt_atgrade','lrt_gradesep']),
                  'mrt_atgrade': rt_isochrones_latlon['rt_mode'] == 'mrt_atgrade',
                  'mrt_gradesep': rt_isochrones_latlon['rt_mode'] == 'mrt_gradesep',
                  'mrt': rt_isochrones_latlon['rt_mode'].isin(['mrt_atgrade','mrt_gradesep']),
                  'all_atgrade': rt_isochrones_latlon['rt_mode'].isin(['brt_atgrade','lrt_atgrade','mrt_atgrade']),
                  'all_gradesep': rt_isochrones_latlon['rt_mode'].isin(['brt_gradesep','lrt_gradesep','mrt_gradesep']),
                  'all': rt_isochrones_latlon['rt_mode'] != None,
                  }
              for mode in list(mode_selectors.keys()):
                  mode_selector = mode_selectors[mode]
                  opened_before = rt_isochrones_latlon['year_open'] <= year
                  not_closed = (np.isnan(rt_isochrones_latlon.year_clos) | (rt_isochrones_latlon.year_clos>year))
                  selector = mode_selector & opened_before & not_closed
                  total_isochrone = gpd.GeoDataFrame(
                      geometry=[rt_isochrones_latlon[selector].unary_union],
                      crs=4326)
                  
                  #fix clipping errors in total_isochone
                  if total_isochrone.unary_union and total_isochrone.unary_union.type == 'GeometryCollection':
                      new_geoms = [x for x in total_isochrone.unary_union.geoms if x.type in ['Polygon','MultiPolygon']]
                      total_isochrone = gpd.GeoDataFrame(geometry=new_geoms,crs=4326)
                  
                  total_isochrone.to_file(f'{folder_name}geodata/rapid_transit/{year}/{mode}_isochrones_ll.geojson', driver='GeoJSON')
                  select_stns = gpd.GeoDataFrame(
                      geometry=[rt_stns[selector].unary_union],
                      crs=4326)
                  select_stns.to_file(f'{folder_name}geodata/rapid_transit/{year}/{mode}_stations_ll.geojson', driver='GeoJSON')
                  
              line_mode_selectors = {
                  'brt_atgrade': rt_lines['rt_mode'] == 'brt_atgrade',
                  'brt_gradesep': rt_lines['rt_mode'] == 'brt_gradesep',
                  'brt': rt_lines['rt_mode'].isin(['brt_atgrade','brt_gradesep']),
                  'lrt_atgrade': rt_lines['rt_mode'] == 'lrt_atgrade',
                  'lrt_gradesep': rt_lines['rt_mode'] == 'lrt_gradesep',
                  'lrt': rt_lines['rt_mode'].isin(['lrt_atgrade','lrt_gradesep']),
                  'mrt_atgrade': rt_lines['rt_mode'] == 'mrt_atgrade',
                  'mrt_gradesep': rt_lines['rt_mode'] == 'mrt_gradesep',
                  'mrt': rt_lines['rt_mode'].isin(['mrt_atgrade','mrt_gradesep']),
                  'all_atgrade': rt_lines['rt_mode'].isin(['brt_atgrade','lrt_atgrade','mrt_atgrade']),
                  'all_gradesep': rt_lines['rt_mode'].isin(['brt_gradesep','lrt_gradesep','mrt_gradesep']),
                  'all': rt_lines['rt_mode'] != None,
                  }
              for mode in list(line_mode_selectors.keys()):
                  mode_selector = line_mode_selectors[mode]
                  opened_before = rt_lines['year_open'] <= year
                  not_closed = (np.isnan(rt_lines.year_clos) | (rt_lines.year_clos>year))
                  selector = mode_selector & opened_before & not_closed
                  select_lines = gpd.GeoDataFrame(
                      geometry=[rt_lines[selector].unary_union],
                      crs=4326)
                  select_lines.to_file(f'{folder_name}geodata/rapid_transit/{year}/{mode}_lines_ll.geojson', driver='GeoJSON')
             
      print(debugcounter); debugcounter+=1 
      if 'pnst' in to_test:
          if not os.path.exists(f'{folder_name}geodata/pnst/'):
              os.mkdir(f'{folder_name}geodata/pnst/')
              
          try:
              protectedbike = gpd.read_file(f"{folder_name}geodata/pnpb/pnpb_latlon_{current_year}.geojson")
              if protectedbike.unary_union is None:
                  protectedbike = gpd.GeoDataFrame(geometry = [], crs=4326)
          except:
              protectedbike = gpd.GeoDataFrame(geometry = [], crs=4326)
          
          try:
              rapidtransport = gpd.read_file(f'{folder_name}geodata/rapid_transit/{current_year}/all_isochrones_ll.geojson')
              if rapidtransport.unary_union is None:
                  rapidtransport = gpd.GeoDataFrame(geometry = [], crs=4326)
          except:
              rapidtransport = gpd.GeoDataFrame(geometry = [], crs=4326)
          
          try:
              frequenttransport = gpd.read_file(f"{folder_name}geodata/pnft/pnft_latlon_{current_year}.geojson")
              if frequenttransport.unary_union is None:
                  frequenttransport = gpd.GeoDataFrame(geometry = [], crs=4326)
          except:
              frequenttransport = gpd.GeoDataFrame(geometry = [], crs=4326)
          
          if protectedbike.unary_union is None: 
              #no bike
              transport_and_bike_latlon = gpd.GeoDataFrame(geometry = [], crs=4326)
          elif rapidtransport.unary_union is None and frequenttransport.unary_union is None:
              #bike, but neither kind of transit
              transport_and_bike_latlon = gpd.GeoDataFrame(geometry = [], crs=4326)
          elif rapidtransport.unary_union is None:
              #only frequent, no rapid
              transport_and_bike_latlon = frequenttransport.overlay(protectedbike, how="intersection")
          elif frequenttransport.unary_union is None:
              #only rapid, no frequent
              try:
                  transport_and_bike_latlon = rapidtransport.overlay(protectedbike, how="intersection")
              except TypeError:
                  newrt = gpd.GeoDataFrame(geometry = [x for x in rapidtransport.geometry[0].geoms if x.type in ['Polygon','MultiPolygon']], crs=rapidtransport.crs)
                  transport_and_bike_latlon = gpd.GeoDataFrame(geometry=newrt, crs=4326).overlay(protectedbike, how="intersection")
          else:
              #all of the above
              rapid_or_frequent = rapidtransport.overlay(frequenttransport, how="union")
              transport_and_bike_latlon = rapid_or_frequent.overlay(protectedbike, how="intersection")
          
          
          transport_and_bike_utm = transport_and_bike_latlon.to_crs(utm_crs)
          if transport_and_bike_utm.geometry.area.sum() != 0:
              transport_and_bike_utm = gpd.overlay(transport_and_bike_utm ,boundaries_utm, how='intersection')
              new_geoms = transport_and_bike_utm.geometry.simplify(services_simplification).make_valid().unary_union
              
              if new_geoms.type == 'GeometryCollection':
                  select_geoms = [x for x in new_geoms.geoms if x.type in ['Polygon','MultiPolygon']]
                  transport_and_bike_latlon= gpd.GeoDataFrame(geometry=select_geoms,crs=4326)
              
              transport_and_bike_latlon = transport_and_bike_utm.to_crs(epsg=4326)
              transport_and_bike_latlon.to_file(f"{folder_name}geodata/pnst/pnst_latlon_{current_year}.geojson", driver='GeoJSON') 
         
              
  print(debugcounter); debugcounter+=1
  
  if 'blocks' in to_test:
      #export all blocks
      blocks_utm = gpd.GeoDataFrame(geometry=[block[0] for block in outblocks], crs=utm_crs)
      blocks_utm['area_utm'] = [block[1] for block in outblocks]
      blocks_utm['perim'] = [block[2] for block in outblocks]
      blocks_utm['oblongness'] = [block[3] for block in outblocks]
      blocks_utm['density'] = [1000000/block[1] for block in outblocks]
      
      filtered_blocks_utm = blocks_utm[
          (blocks_utm.oblongness < 50) &
          (blocks_utm.area_utm > 1000) &
          (blocks_utm.area_utm < 1000000)]
      
      #TODO -- determine whether to normal simplify or toposimplify, and how much
      filtered_blocks_utm.geometry = filtered_blocks_utm.geometry.simplify(10)
      blocks_latlon = filtered_blocks_utm.to_crs(epsg=4326)
      #blocks_topo = topojson.Topology(blocks_latlon, prequantize=True)
      #blocks_latlon = blocks_topo.toposimplify(blocks_simplification).to_gdf()
      
      blocks_latlon.to_file(f"{folder_name}geodata/blocks/blocks_latlon_{current_year}.geojson", driver='GeoJSON')
      
      block_patches_latlon, block_unbuffered_patches_latlon = make_patches(
          boundaries_latlon, 
          utm_crs, 
          patch_length=block_patch_length
          )
      block_unbuf_patches_utm = block_unbuffered_patches_latlon.to_crs(utm_crs)
      #export            
      patch_densities = block_unbuffered_patches_latlon
      for patch_idx  in list(patch_densities.index):
          try:
              patch_densities.loc[patch_idx,'block_count'] = blocks_latlon.intersects(patch_densities.loc[patch_idx,'geometry'].centroid).value_counts()[True]
          except KeyError:
              patch_densities.loc[patch_idx,'block_count'] = 0 
      patch_densities_utm = patch_densities.to_crs(utm_crs)
      patch_densities_utm['density'] = patch_densities_utm.block_count / (patch_densities_utm.area / 1000000)
      patch_densities_latlon = patch_densities_utm.to_crs(epsg=4326)
      patch_densities_latlon.to_file(f"{folder_name}geodata/blocks/block_densities_latlon_{current_year}.geojson", driver='GeoJSON')

      
                  
