import warnings
import datetime
import os
import os.path
import numpy
import math
import statistics
import rasterstats
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

import isochrones
import get_service_locations
import gtfs_parser
#import prep_bike_osm
#import prep_pop_ghsl
#import summarize_ttm

import pdb


def make_patches(boundaries_latlon, crs_utm, patch_length = 10000, buffer = 500): #patch_length and buffer in m
    ''' 
    'tile' the boundaries of a city into patches, like a patchwork quilt, including a buffer
    buffer should be half the length of the longest PNX distance
    '''
    bounds_utm_poly = boundaries_latlon.to_crs(crs_utm).geometry.unary_union
    bbox_utm = bounds_utm_poly.bounds
    
    height_m = abs(bbox_utm[3]-bbox_utm[1])
    width_m = abs(bbox_utm[2]-bbox_utm[0])
    
    n_hslicers = math.floor(height_m / patch_length)
    n_vslicers = math.floor(width_m / patch_length)
    
    hslicers = []
    vslicers = []
    
    for i in range(1,n_hslicers+1):
        h_increment = height_m/(n_hslicers+1)
        lat = bbox_utm[1]+(i*h_increment)
        slicer = shapely.geometry.LineString([(bbox_utm[0],lat),(bbox_utm[2],lat)])
        hslicers.append(slicer)
        
    for i in range(1,n_vslicers+1):
        v_increment = width_m/(n_vslicers+1)
        lon = bbox_utm[0]+(i*v_increment)
        slicer = shapely.geometry.LineString([(lon,bbox_utm[1]),(lon, bbox_utm[3])])
        vslicers.append(slicer)
    
    
    patches = [bounds_utm_poly] #a list of polygons
    
    for slicer in hslicers+vslicers:
        newpatches = []
        for patch in patches:
            newpatches += list(shapely.ops.split(patch, slicer).geoms)
        patches = newpatches
    
    buffered_patches = []
    for patch in patches:
        buffered_patch = patch.buffer(buffer).intersection(bounds_utm_poly)
        buffered_patches.append(buffered_patch)
    patches_utm = gpd.GeoDataFrame(geometry=buffered_patches, crs=crs_utm)
    patches_latlon = patches_utm.to_crs(4326)
    patches_latlon['patch_id'] = patches_latlon.index
    unbuffered_patches_utm = gpd.GeoDataFrame(geometry=patches, crs=crs_utm)
    unbuffered_patches_latlon = unbuffered_patches_utm.to_crs(4326)
    unbuffered_patches_latlon['patch_id'] = patches_latlon.index
    
    print(f"cut {len(patches)} patches")
    
    return patches_latlon, unbuffered_patches_latlon

def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    #TODO: could I make this more efficient by taking the intersections of circles
    # with a radius of distance, distance*2, etc?
    #assuming that in most cases the line will be relatively straight.
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                [LineString(coords[:i+1]),
                LineString(coords[i:])],
                pd]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                [LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])],
                cp]

def add_intermediate_nodes(G, maxdist, exceptionsize=1000, debug=False, folder_name=''): #TODO CUT?
    start_time = datetime.datetime.now()
    print(f'adding nodes to a max edge length of {maxdist}')
    nodes, edges = ox.graph_to_gdfs(G)
    newnode_baseid = nodes.index.max()+1
    nodes_added = 0
    # new_nodes = gpd.GeoDataFrame(crs=nodes.crs, geometry=[])
    # new_edges = gpd.GeoDataFrame(crs=nodes.crs, geometry=[])
    # new_edges.index.astype(edges.index.dtype, copy=False)
    too_long = edges[(edges.geometry.length>maxdist) &
                     (edges.geometry.length<exceptionsize)] # we leave out ways > 1km
    for to_shorten in tqdm(too_long.iterrows(), total = len(too_long)):
        to_shorten = to_shorten[1]
        old_leg = to_shorten.geometry
        
        while old_leg.length > maxdist:
            new_legs, new_node = cut(old_leg, maxdist/2)
        
            new_node_id = newnode_baseid + nodes_added
            nodes_added += 1
            nodes.loc[new_node_id, 'y'] = new_node.y
            nodes.loc[new_node_id, 'x'] = new_node.x
            nodes.loc[new_node_id,'osmid'] = f'NEW_NODE_{new_node_id}'
            
            first_leg_id = (to_shorten.name[0], new_node_id, to_shorten.name[2])
            edges.loc[first_leg_id,'length'] = new_legs[0].length
            edges.loc[first_leg_id,'geometry'] = new_legs[0]
            edges.loc[first_leg_id,'osmid'] = f'NEW_EDGE_{new_node_id}'
            
            old_leg = new_legs[1]
            
        second_leg_id = (new_node_id, to_shorten.name[1], to_shorten.name[2])
        edges.loc[second_leg_id,'geometry'] = new_legs[1]
        edges.loc[second_leg_id,'length'] = new_legs[1].length
        edges.loc[second_leg_id,'osmid'] = f'NEW_EDGE_{new_node_id+1}'
        
        edges = edges.drop(to_shorten.name)
        
    
    end_time = datetime.datetime.now()
    print(f"added {nodes_added} nodes in {end_time - start_time}")
    return ox.graph_from_gdfs(nodes, edges)

        

def weighted_pop_density(array):
    total = 0
    for cell in array.compressed():
        total += cell**2
    return total / numpy.sum(array)

def spatial_analysis(boundaries, 
                      id_code,
                      name,
                      folder_name='', 
                      buffer_dist=100,#buffer out from nodes, m
                      headway_threshold=10,#min
                      to_test = [
                           'healthcare',
                           'schools',
                           'h+s',
                           #'libraries',
                           'bikeshare',
                           'carfree',
                           'blocks',
                           'density',
                           'pnft',
                           'pnrt',
                           'pnpb', #protected bikeways
                           'pnab', #all bikeways
                           'highways',
                           #'access',
                           #'transport_performance',
                           #'special',
                           ],
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
                      years = range(1975,2031), #for PNRT and pop_dens. remember range(1,3) = [1,2]
                      current_year = 2022,
                      overpass = False,
                      patch_length = 8000, #m
                      block_patch_length = 2000, #m
                      boundary_buffer = 500, #m
                      blocks_simplification = 0.0001, #topo-simplification
                      services_simplification = 10, #m TODO replace with gpd simplification?
                      access_resolution = 1000, #m
                      transport_performance_speeds = { #make sure these align with the r script
                          'walk': 3.6, #km/h
                          'bike_lts1': 12,
                          'bike_lts2': 12,
                          #'transit': 12, #more arbitrary
                          },
                      transport_performance_times = [30,45,60], #mins
                      gtfs_files = [],
                      ghsl_projection = 'mw', #TODO make this do something
                      ghsl_resolution = '1000',
                      #max_edge_length = , #should not be bigger than 2xbuffer_dist
                      debug = True,
                      ):   
    
    dt = datetime.datetime.now()
    
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    ox.utils.config(log_console = False)
    
    useful_tags = ox.settings.useful_tags_way + ['cycleway', 'cycleway:left', 'cycleway:right', 'cycleway:both', 'bicycle']
    ox.config(use_cache=True, log_console=True, useful_tags_way=useful_tags)
    walk_filter = ('["area"!~"yes"]["highway"]["highway"!~"link|motor'
                   '|proposed|construction|abandoned'
                   '|platform|raceway"]'
                   '["service"!~"parking_aisle|driveway"]'
                   '["foot"!~"no"]["service"!~"private"]'
                   '{}').format(ox.settings.default_access) #used later
    
    if folder_name != '' and not folder_name[-1:] == '/':
        folder_name += '/'
    if 'h+s' in to_test:
        for part in ['healthcare', 'schools']:
            if not part in to_test:
                to_test.append(part)
                    
    
    boundaries_latlon = gpd.GeoDataFrame(geometry = [boundaries])
    boundaries_latlon.crs = {'init':'epsg:4326'}
    longitude = round(numpy.mean(boundaries_latlon.geometry.centroid.x),10)
    utm_zone = int(math.floor((longitude + 180) / 6) + 1)
    utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
    boundaries_utm = boundaries_latlon.to_crs(utm_crs)
    boundaries_utm.geometry = boundaries_utm.geometry.buffer(boundary_buffer)
    boundaries_latlon = boundaries_utm.to_crs(epsg=4326)
    boundaries_mw = boundaries_utm.to_crs("ESRI:54009")
    boundaries = boundaries_latlon.geometry.unary_union
    
    print('Calculating geodata for Pedestrians First indicators in',name)
    print('Measuring',str(to_test))
    
    for year in years:
        if year % 5 == 0:
            if year < current_year:
                in_file = f'input_data/ghsl/GHS_POP_E{year}_GLOBE_R2022A_54009_{ghsl_resolution}/GHS_POP_E{year}_GLOBE_R2022A_54009_{ghsl_resolution}_V1_0.tif'
            else:
                in_file = f'input_data/ghsl/GHS_POP_P{year}_GLOBE_R2022A_54009_{ghsl_resolution}/GHS_POP_P{year}_GLOBE_R2022A_54009_{ghsl_resolution}_V1_0.tif'
            with rasterio.open(in_file) as dataset:
                out_image, out_transform = rasterio.mask.mask(dataset, [boundaries_mw.unary_union], crop=True)
                out_meta = dataset.meta
                out_meta.update({"driver": "GTiff",
                      "height": out_image.shape[1],
                      "width": out_image.shape[2],
                      "transform": out_transform})
            if not os.path.exists(f"{folder_name}geodata/population"):
                os.mkdir(f"{folder_name}geodata/population")
            with rasterio.open(f"{folder_name}geodata/population/pop_{year}.tif", "w", **out_meta) as dest:
                dest.write(out_image)
    
    testing_services = []
    for service in ['healthcare', 'schools', 'libraries', 'bikeshare']:
        if service in to_test:
            testing_services.append(service)
            
    service_point_locations={}
    if len(testing_services) > 0:
        if overpass:
            raise RuntimeError
            #haven't set this up yet, because I don't have a query for bikeshare or carfree
        else:
            handler = get_service_locations.ServiceHandler()
            handler.apply_file(folder_name+'temp/city.o5m', locations=True)
            for service in testing_services:
                coords = handler.locationlist[service]
                service_point_locations[service] = gpd.GeoDataFrame(
                    geometry = [Point(coord) for coord in coords],
                    crs=4326)
            citywide_carfree = handler.carfreelist
                
    if 'pnab' in to_test or 'pnpb' in to_test:
        testing_services.append('pnpb')
        testing_services.append('pnab')
        
    if 'highways' in to_test:
        all_highway_polys = []
            
    #should I finish this?
    # if 'highways' in to_test:     
    #     highway_filter = '["area"!~"yes"]["highway"~"motorway|trunk"]'
    #     G_hwys = ox.graph_from_polygon(boundaries,
    #                             custom_filter=highway_filter,
    #                             retain_all=True)
    #     G_hwys = ox.project_graph(G_hwys)
    #     edges_hwys = ox.utils_graph.graph_to_gdfs(G_hwys, nodes=False)
    #     edges_polyline = edges_hwys.geometry.unary_union
    #     increment_dist = 25
    #     distances = np.arange(0, edges_polyline.length, increment_dist)
    #     points = [edges_polyline.interpolate(distance) for distance in distances] + [edges_polyline.boundary[1]]
    #     multipoint = unary_union(points)
    #     #need to switch back to latlon and put in all_coords['highway']
        
    if 'special' in to_test:
        special = gpd.read_file(f'{folder_name}special.geojson')
        if special.geometry[0].type == 'MultiLineString':
            special.crs = 4326
            special_utm = special.to_crs(utm_crs)
            points=[]
            for line in special_utm.geometry.unary_union:
                line = shapely.ops.transform(lambda x, y, z=None: (x, y), line)
                while line.length > 50:
                    cutlines, point = cut(line, 50)
                    points.append(point)
                    line=cutlines[1]
            points_utm = gpd.GeoDataFrame(geometry=points, crs=utm_crs)
            points_latlon = points_utm.to_crs(4326)
            service_point_locations['special'] = points_latlon
        
        elif special.geometry[0].type == 'LineString':
            special.crs = 4326
            special_utm = special.to_crs(utm_crs)
            points=[]
            for line in special_utm.geometry:
                line = shapely.ops.transform(lambda x, y, z=None: (x, y), line)
                while line.length > 50:
                    cutlines, point = cut(line, 50)
                    points.append(point)
                    line=cutlines[1]
            points_utm = gpd.GeoDataFrame(geometry=points, crs=utm_crs)
            points_latlon = points_utm.to_crs(4326)
            service_point_locations['special'] = points_latlon
        elif special.geometry[0].type == 'Point':
            service_point_locations['special'] = special
        else:
            import pdb; pdb.set_trace()
            raise RuntimeError
        testing_services.append('special')

    if 'pnft' in to_test:
        testing_services.append('pnft')
        freq_stops = gtfs_parser.get_frequent_stops(
            boundaries, 
            folder_name, 
            headway_threshold)
        service_point_locations['pnft'] = freq_stops
            
    if 'pnrt' in to_test:
        mode_classifications = {
            'Bus Rapid Transit':'brt',
            'Light Rail': 'lrt',
            'Light Metro': 'lrt',
            'Heavy Rail': 'mrt',
            }
        rt_lines = gpd.read_file('input_data/transit_explorer/shapefile/lines.shp')
        rt_stns = gpd.read_file('input_data/transit_explorer/shapefile/stations.shp')
        rt_lines = rt_lines.overlay(boundaries_latlon, how='intersection')
        rt_stns = rt_stns.overlay(boundaries_latlon, how='intersection')
        rt_lines = rt_lines[rt_lines['mode'].isin(mode_classifications.keys())]
        rt_stns = rt_stns[rt_stns['mode'].isin(mode_classifications.keys())]
        for idx in rt_lines.index:
            rt_lines.loc[idx,'rt_mode'] = mode_classifications[rt_lines.loc[idx,'mode']]
        for idx in rt_stns.index:
            rt_stns.loc[idx,'rt_mode'] = mode_classifications[rt_stns.loc[idx,'mode']]
        rt_isochrones = rt_stns.copy()
        rt_stns_utm = rt_stns.to_crs(utm_crs)
        rt_isochrones_utm = rt_isochrones.to_crs(utm_crs)
    
    # if 'access' in to_test or 'transport_performance' in to_test:
        
        # for file in [folder_name+'temp/access/grid_pop.geojson',
        #              folder_name+'temp/access/city_ltstagged.pbf']:
        #     if os.path.exists(file):
        #         os.remove(file)
        
        # #prep osm (add LTS values)
        # original_filename = folder_name+"temp/city.pbf"
        # prep_bike_osm.add_lts_tags(original_filename,folder_name+"temp/access/city_ltstagged.pbf")
        
        # #prep pop
        # prep_pop_ghsl.setup_grid(
        #     utm_crs, 
        #     boundaries_utm.unary_union, 
        #     access_resolution, 
        #     folder_name+"geodata/pop_dens.tif", 
        #     adjust_pop = True,
        #     save_loc = (folder_name+'temp/access/pop_points.csv',
        #                 folder_name+'temp/access/grid_pop.geojson')
        #     )
        # grid_pop = gpd.read_file(folder_name+'temp/access/grid_pop.geojson')
        
        # #cp over script
        # shutil.copy('access/two_step_access/calcttm_simple.r', folder_name+'temp/access/calcttm_simple.r')
        
        # #run script
        # command = f"Rscript calcttm_simple.r {folder_name}temp/access/"
        # subprocess.check_call(command.split(' '))
        
        # #
        # all_modes = ['walk','car','pnft','bike_lts1','bike_lts2','bike_lts4']
        
        # mode_ttms = {}
        # actual_modes = []
        # for mode in all_modes:
        #     if os.path.exists(f'{folder_name}temp/access/{mode}_wide.csv'):
        #         mode_ttms[mode] = summarize_ttm.load_wide_ttm(f'{folder_name}temp/access/{mode}_wide.csv')
        #         actual_modes.append(mode)
        
        
        # total_value, cxn_df_latlon = summarize_ttm.evaluation(grid_pop, 
        #                                                mode_ttms, 
        #                                                cumulative_time_lims = transport_performance_times, 
        #                                                ec_modes = [], 
        #                                                ec_ttms = None)
        
        # if 'transport_performance' in to_test:
        #     cxn_df_utm = ox.project_gdf(cxn_df_latlon)
        #     for mode in transport_performance_speeds.keys():
        #         for lim in transport_performance_times:
        #             tp_buffer_dist = transport_performance_speeds[mode] * (1000/60) * lim
        #             cxn_df_utm[f'buffer_{mode}_{lim}'] = cxn_df_utm.centroid.buffer(tp_buffer_dist)
        #     for point_id in cxn_df_utm.index:
        #         for mode in transport_performance_speeds.keys():
        #             for lim in transport_performance_times:
        #                 actual_sum = cxn_df_utm.loc[point_id, f'cumsum_{mode}_{lim}']
        #                 circle = cxn_df_utm.loc[point_id, f'buffer_{mode}_{lim}']
        #                 circ_gdf = gpd.GeoDataFrame(geometry = [circle], crs = cxn_df_utm.crs)
        #                 potential_sum = cxn_df_utm.overlay(circ_gdf, how='intersection').population.sum()
        #                 cxn_df_utm.loc[point_id, f'potentialsum_{mode}_{lim}'] = potential_sum
        #                 cxn_df_utm.loc[point_id, f'performance_{mode}_{lim}'] = actual_sum / potential_sum
        #     # before saving the file convert geometry to wkt
        #     cxn_df_latlon = cxn_df_utm.to_crs(4326)
        #     for mode in transport_performance_speeds.keys():
        #         for lim in transport_performance_times:
        #             cxn_df_latlon[f'buffer_{mode}_{lim}'] = cxn_df_utm[f'buffer_{mode}_{lim}'].to_crs(4326).to_wkt()
            
        # #save connections_df
        # cxn_df_latlon.to_file(f'{folder_name}geodata/connections.gpkg', driver='GPKG')
    
    
    quilt_isochrone_polys = {}
    for service in to_test:
        quilt_isochrone_polys[service] = False
        if 'pnab' in to_test or 'pnpb' in to_test:
            quilt_isochrone_polys['pnab'] = False
            quilt_isochrone_polys['pnpb'] = False
            quilt_allbike = gpd.GeoDataFrame(geometry=[], crs=utm_crs)
            quilt_protectedbike = gpd.GeoDataFrame(geometry=[], crs=utm_crs)
    
    
    patches, unbuffered_patches = make_patches(boundaries_latlon, utm_crs, patch_length=patch_length)
    if debug:
        patches.to_file(folder_name+'debug/patches.geojson', driver='GeoJSON')
    
    patch_times = []
    if len(to_test) > 0 and to_test != ["blocks"]:
        for p_idx in patches.index:
            patch = patches.loc[p_idx, 'geometry']
            unbuffered_patch = unbuffered_patches.loc[p_idx, 'geometry']
            try:
                patch_start = datetime.datetime.now()
                for cleanupfile in ['pathbounds.geojson',
                                    'pathbounds.poly',
                                    'patch.osm', 
                                    'allhwyspatch.osm']:
                    if os.path.exists('temp/'+cleanupfile):
                        os.remove('temp/'+cleanupfile)
                
                patchgdf = gpd.GeoDataFrame(geometry=[patch], crs=4326)
                
                patchgdf.to_file('temp/patchbounds.geojson', driver='GeoJSON')
                subprocess.run('python ogr2poly/ogr2poly.py temp/patchbounds.geojson > temp/patchbounds.poly', shell=True, check=True)
                
                if overpass:
                    G = ox.graph_from_polygon(patch, 
                                              custom_filter=walk_filter, 
                                              simplify=True, 
                                              retain_all=True)
                else:
                    try:
                        subprocess.check_call(['osmconvert',
                                               str(folder_name)+'temp/cityhighways.o5m',
                                               "-B=temp/patchbounds.poly",
                                               #'--complete-ways',  #was commented
                                               '--drop-broken-refs',  #was uncommented
                                               '-o=temp/patch_allroads.osm'])
                        G_allroads = ox.graph_from_xml('temp/patch_allroads.osm', simplify=True, retain_all=True)
                        os.remove('temp/patch_allroads.osm')
                        
                        # TESTING using the same graph for everything. 
                        # A little more efficient, avoids some problems with cycleways
                        # But implies that people can walk along motorways.
                        # Hopefully not a problem since those are usually divided.
                        # If I wanted to get really fancy I could take them out 
                        # entirely, before doing isochrones, I guess
                        G = G_allroads
                        
                        # subprocess.check_call(['osmconvert',
                        #                        str(folder_name)+'temp/citywalk.o5m',
                        #                        "-B=temp/patchbounds.poly",
                        #                        #'--complete-ways',  #was commented
                        #                        '--drop-broken-refs',  #was uncommented
                        #                        '-o=temp/patch.osm'])
                        # G = ox.graph_from_xml('temp/patch.osm', simplify=True, retain_all=True)
                        # os.remove('temp/patch.osm')
                    except TypeError: #something to do with clipping, seems to happen once in a while
                        #pdb.set_trace()
                        #this is a very stupid band-aid, but it works for now
                        #TODO either figure this out or make the patch more efficient somehow? idk
                        #I think right now it's getting EVERYTHING, not just highways. Except not highways that are abandoned :)
                        print ('TYPEERROR FROM CLIPPING PATCH', p_idx)
                        with open(str(folder_name)+"patcherrorlog.txt", "a") as patcherrorlog:
                            patcherrorlog.write('TYPEERROR FROM CLIPPING PATCH '+str(p_idx))
                        G = ox.graph_from_polygon(patch, 
                                              #custom_filter=walk_filter, 
                                              simplify=True, 
                                              retain_all=True)
                os.remove('temp/patchbounds.poly')
                        
                G.remove_nodes_from(list(nx.isolates(G)))
                
                if len(G.edges) > 0 and len(G.nodes) > 0:
                    G = ox.project_graph(G, to_crs=utm_crs)
                    G_allroads = ox.project_graph(G_allroads, to_crs=utm_crs)
                    
                    center_nodes = {}
                    if 'pnab' in to_test or 'pnpb' in to_test:
                        ways_gdf = ox.graph_to_gdfs(G_allroads, nodes=False)
                        
                        for col in ['highway','cycleway','bicycle','cycleway:left','cycleway:right','cycleway:both']:
                            if not col in ways_gdf.columns:
                                ways_gdf[col] = ''
                                
                        tagged_cycleways = ways_gdf[(ways_gdf['highway'] == 'cycleway')]
                        cycle_paths = ways_gdf[(ways_gdf['highway'] == 'path') & (ways_gdf['bicycle'] == 'designated')]
                        on_street_tracks = ways_gdf[(ways_gdf['cycleway:right'] == 'track') |
                                                       (ways_gdf['cycleway:left'] == 'track') |
                                                       (ways_gdf['cycleway:both'] == 'track') |
                                                       (ways_gdf['cycleway'] == 'track')]
                        on_street_lanes = ways_gdf[(ways_gdf['cycleway:right'] == 'lane') |
                                                       (ways_gdf['cycleway:left'] == 'lane') |
                                                       (ways_gdf['cycleway:both'] == 'lane') |
                                                       (ways_gdf['cycleway'] == 'lane')]
                        
                        total_protectedbike = pd.concat([tagged_cycleways, cycle_paths, on_street_tracks])
                        total_allbike = pd.concat([total_protectedbike, on_street_lanes])
                        
                        quilt_allbike = pd.concat([quilt_allbike, total_allbike])
                        quilt_protectedbike = pd.concat([quilt_protectedbike, total_protectedbike])
                        
                        center_nodes['pnpb'] = set()
                        center_nodes['pnab'] = set()
                        for edge in total_protectedbike.index:
                            center_nodes['pnpb'].add(edge[0])
                            center_nodes['pnpb'].add(edge[1])
                        for edge in total_allbike.index:
                            center_nodes['pnab'].add(edge[0])
                            center_nodes['pnab'].add(edge[1])
                    
                    
                    if 'highways' in to_test:
                        highway_polys = get_service_locations.get_highways(G_allroads)
                        if highway_polys is not None:
                            all_highway_polys.append(highway_polys)
                        
                    # if debug:
                    #     for node, data in G.nodes(data=True):
                    #         if 'osmid' in data:
                    #             data['osmid_original'] = data.pop('osmid')
                    #     ox.io.save_graph_geopackage(G, f'{folder_name}debug/{p_idx}_graph.gpkg')
                        
                    for service in service_point_locations.keys():
                        if service in ['healthcare','schools','libraries','bikeshare','pnft','special']:
                            patch_points = service_point_locations[service].intersection(patch)
                            patch_points.crs=4326
                            patch_points_utm = patch_points.to_crs(utm_crs)
                            patch_points_utm = patch_points_utm[(~patch_points_utm.is_empty)
                                                                &(patch_points_utm.geometry!=None)]
                            if len(patch_points) > 0:
                                center_nodes[service]  = ox.distance.nearest_nodes(
                                    G, 
                                    patch_points_utm.geometry.x, 
                                    patch_points_utm.geometry.y,
                                    return_dist=False)
                    # if debug:
                    #     with open(f'{folder_name}debug/{p_idx}_centernodes.json', "w") as i :
                    #         json.dump(center_nodes, i)
                    
                    
                    isochrone_polys = {}
                    failures = {}
                    for service in to_test:
                        failures[service] = 0
                        if 'pnab' in to_test or 'pnpb' in to_test:
                            failures['pnab']=0
                            failures['pnpb']=0
                    
                    # Get polygons
                    for service in testing_services:
                        if service in center_nodes.keys():
                            print(f'getting polygons for {service}, {len(center_nodes[service])} center_nodes')
                            isochrone_polys[service] = isochrones.proper_iso_polys(
                                G, 
                                center_nodes[service],
                                distance=distances[service],                                           
                                buffer=buffer_dist, 
                                infill=2500)       
                        
                    for service in isochrone_polys.keys():
                        if service not in quilt_isochrone_polys.keys() or not quilt_isochrone_polys[service]:
                            quilt_isochrone_polys[service] = isochrone_polys[service]
                        elif isochrone_polys[service]:
                            quilt_isochrone_polys[service] = shapely.ops.unary_union([quilt_isochrone_polys[service],isochrone_polys[service]])
                        if debug == True:
                            if not os.path.exists(f'{folder_name}debug/patch{p_idx}'):
                                os.makedirs(f'{folder_name}debug/patch{p_idx}')
                            gpd.GeoSeries([isochrone_polys[service]], crs=utm_crs).to_crs(4326).to_file(f'{folder_name}debug/patch{p_idx}/{service}_latlon.geojson', driver='GeoJSON')
                    
                    #get polygons for rapid transit
                    if 'pnrt' in to_test:
                        for stn_idx in rt_stns.index:
                            if unbuffered_patch.contains(rt_stns.loc[stn_idx, 'geometry']):
                                stn_utm = rt_stns_utm.loc[stn_idx, 'geometry']
                                center_node = ox.distance.nearest_nodes(
                                    G, 
                                    stn_utm.x, 
                                    stn_utm.y,
                                    return_dist=False)
                                iso_poly = isochrones.proper_iso_polys(
                                    G, 
                                    [center_node],
                                    distance=distances['pnrt'],                                           
                                    buffer=buffer_dist, 
                                    infill=2500)
                                rt_isochrones_utm.loc[stn_idx,'geometry'] = iso_poly
                    
                    #import pdb;pdb.set_trace()
                    patch_time = datetime.datetime.now() - patch_start
                        
                    print(f"finished patch #{p_idx+1} out of {len(patches)} in {patch_time}")
                    patch_times.append(patch_time)
            except ox._errors.EmptyOverpassResponse:
                print('EmptyOverpassResponse')
                pass
            #except:
            #    print("GOT SOME ERROR FOR PATCH", p_idx)
            # except ValueError:
            #     print('ValueError')
            #     now = str(datetime.datetime.now())
            #     with open('error'+now+'.txt','w') as errout:
            #         traceback.print_exc(limit=3,file=errout)
            #     print('saved to error'+now+'.txt') 
    
        if debug:
            patches.times = patch_times
            patches.to_file(folder_name+'debug/patches_with_times.geojson', driver='GeoJSON')
            pd.DataFrame({'failures':failures}).to_csv(folder_name+'debug/failures.csv')
    
    #start saving files
    for service in testing_services:
        if quilt_isochrone_polys[service]:
            service_utm = gpd.GeoDataFrame(geometry = [quilt_isochrone_polys[service]],
                                           crs=utm_crs)
            service_utm.geometry = service_utm.geometry.simplify(services_simplification)
            service_utm = gpd.overlay(service_utm ,boundaries_utm, how='intersection')
            service_latlon = service_utm.to_crs(epsg=4326)
            #service_utm.to_file(folder_name+service+'utm'+'.geojson', driver='GeoJSON')
            service_latlon.to_file(folder_name+'geodata/'+service+'latlon'+'.geojson', driver='GeoJSON')
        if service in service_point_locations.keys():
            service_point_locations[service].to_file(folder_name+'geodata/'+service+'_points_latlon'+'.geojson', driver='GeoJSON')
            
    if 'pnpb' in to_test or 'pnab' in to_test:
        if not quilt_protectedbike.empty:
            quilt_protectedbike = quilt_protectedbike.to_crs(4326)
            merged_protectedbike = quilt_protectedbike.intersection(boundaries)
            merged_protectedbike = gpd.GeoDataFrame(geometry = [merged_protectedbike.unary_union], crs=4326)
            merged_protectedbike.to_file(folder_name+'geodata/protectedbike_latlon.geojson',driver='GeoJSON')
        if not quilt_allbike.empty:
            quilt_allbike = quilt_allbike.to_crs(4326)
            merged_allbike = quilt_allbike.intersection(boundaries)
            merged_allbike = gpd.GeoDataFrame(geometry = [merged_allbike.unary_union], crs=4326)
            merged_allbike.to_file(folder_name+'geodata/allbike_latlon.geojson',driver='GeoJSON')
    
    if 'highways' in to_test:
        all_hwys_poly = shapely.ops.unary_union(all_highway_polys)
        all_hwys_gdf = gpd.GeoDataFrame(geometry=[all_hwys_poly], crs=utm_crs)
        all_hwys_latlon = all_hwys_gdf.to_crs(4326)
        all_hwys_latlon.to_file(folder_name+'geodata/allhwys_latlon.geojson',driver='GeoJSON')
        buffered_hwys = all_hwys_gdf.buffer(distances['highways'])
        buffered_hwys_latlon = buffered_hwys.to_crs(4326)
        buffered_hwys_latlon.to_file(folder_name+'geodata/buffered_hwys_latlon.geojson',driver='GeoJSON')
        
        
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
        carfree_latlon.to_file(folder_name+'geodata/carfreelatlon.geojson', driver='GeoJSON')
        
    if 'h+s' in to_test:
        if quilt_isochrone_polys['healthcare'] and quilt_isochrone_polys['schools']:
            service = 'h+s'
            intersect = shapely.ops.unary_union([quilt_isochrone_polys['healthcare'].intersection(quilt_isochrone_polys['schools'])])
            if type(intersect) == shapely.geometry.collection.GeometryCollection:
                intersect = [obj for obj in intersect if type(obj) == shapely.geometry.polygon.Polygon]
                intersect = shapely.geometry.MultiPolygon(intersect)
            hs_utm = gpd.GeoDataFrame(geometry = [intersect], crs=utm_crs)
            if hs_utm.geometry.area.sum() != 0:
                hs_utm = gpd.overlay(hs_utm ,boundaries_utm, how='intersection')
                hs_utm.geometry = hs_utm.geometry.simplify(services_simplification)
                hs_latlon = hs_utm.to_crs(epsg=4326)
                hs_latlon.to_file(folder_name+'geodata/'+service+'latlon.geojson', driver='GeoJSON')
    
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
                    'brt': rt_isochrones_latlon['rt_mode'] == 'brt',
                    'lrt': rt_isochrones_latlon['rt_mode'] == 'lrt',
                    'mrt': rt_isochrones_latlon['rt_mode'] == 'mrt',
                    'all': rt_isochrones_latlon['rt_mode'].isin(['brt','lrt','mrt']),
                    }
                for mode in ['brt','lrt','mrt','all']:
                    mode_selector = mode_selectors[mode]
                    opened_before = rt_isochrones_latlon['year_open'] <= year
                    not_closed = (np.isnan(rt_isochrones_latlon.year_close) | (rt_isochrones_latlon.year_close>year))
                    selector = mode_selector & opened_before & not_closed
                    total_isochrone = gpd.GeoDataFrame(
                        geometry=[rt_isochrones_latlon[selector].unary_union],
                        crs=4326)
                    total_isochrone.to_file(f'{folder_name}geodata/rapid_transit/{year}/{mode}_isochrones_ll.geojson', driver='GeoJSON')
                    select_stns = gpd.GeoDataFrame(
                        geometry=[rt_stns[selector].unary_union],
                        crs=4326)
                    select_stns.to_file(f'{folder_name}geodata/rapid_transit/{year}/{mode}_stations_ll.geojson', driver='GeoJSON')
                line_mode_selectors = {
                    'brt': rt_lines['rt_mode'] == 'brt',
                    'lrt': rt_lines['rt_mode'] == 'lrt',
                    'mrt': rt_lines['rt_mode'] == 'mrt',
                    'all': rt_lines['rt_mode'].isin(['brt','lrt','mrt']),
                    }
                for mode in ['brt','lrt','mrt','all']:
                    mode_selector = mode_selectors[mode]
                    opened_before = rt_lines['year_open'] <= year
                    not_closed = (np.isnan(rt_lines.year_close) | (rt_lines.year_close>year))
                    selector = mode_selector & opened_before & not_closed
                    select_lines = gpd.GeoDataFrame(
                        geometry=[rt_lines[selector].unary_union],
                        crs=4326)
                    select_lines.to_file(f'{folder_name}geodata/rapid_transit/{year}/{mode}_lines_ll.geojson', driver='GeoJSON')
               
    if 'blocks' in to_test:
        print("getting blocks")
        #TODO -- should probably eventually fix this so that the patches
        #used in getting the blocks are not necessarily the same
        #as the ones used in summarizing them
        patches_latlon, unbuffered_patches_latlon = make_patches(
            boundaries_latlon, 
            utm_crs, 
            patch_length=block_patch_length
            )
        unbuf_patches_utm = unbuffered_patches_latlon.to_crs(utm_crs)
        print ("cut", len(list(patches)),"patches for block size in",name)
        
        outblocks = []
        block_counts = []
        for patch_idx in patches_latlon.index:
            patch_latlon = patches_latlon.loc[patch_idx, 'geometry']
            unbuf_patch_utm = unbuf_patches_utm.loc[patch_idx, 'geometry']
            
            print("patch"+str(patch_idx)+" of "+str(len(patches_latlon)) )
            if overpass:
                G = ox.graph_from_polygon(patch, custom_filter=walk_filter, simplify=True, retain_all=True)
            else:
                boundingarg = '-b='
                boundingarg += str(patch_latlon.bounds[0])+','
                boundingarg += str(patch_latlon.bounds[1])+','
                boundingarg += str(patch_latlon.bounds[2])+','
                boundingarg += str(patch_latlon.bounds[3])
                subprocess.check_call(['osmconvert',
                                       folder_name+'temp/citywalk.o5m',
                                       boundingarg,
                                       #'--complete-ways',
                                       '--drop-broken-refs',
                                       '-o=patch.osm'])
                try:
                    G = ox.graph_from_xml('patch.osm', simplify=True, retain_all=True)
                except:
                    G = False
            
            if G and len(G.edges) > 0:
                G = ox.project_graph(G, to_crs=utm_crs)
                
                streets = ox.utils_graph.graph_to_gdfs(G, nodes = False)
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
                            if block.centroid.within(unbuf_patch_utm):
                                area = round(block.area, 3)
                                perim = round(block.length, 3)
                                lemgth = round((perim * perim) / area, 3)
                                if blocks_simplification:
                                    block = block.simplify(blocks_simplification)
                                all_blocks.append((block, area, perim, lemgth))
                                if (lemgth < 50) and (1000 < area < 1000000):
                                    selected_areas.append(area)
                    outblocks += all_blocks
                    print(f'cut {len(all_blocks)} blocks')
                    block_counts.append(len(all_blocks))
                else:
                    block_counts.append(0)
                    print('not merged!')
            else:
                block_counts.append(0)
        
        #export            
        patch_densities = unbuffered_patches_latlon.copy()
        patch_densities['block_count'] = block_counts
        patch_densities_utm = patch_densities.to_crs(utm_crs)
        patch_densities_utm['density'] = patch_densities_utm.block_count / (patch_densities_utm.area /1000000)
        patch_densities_latlon = patch_densities_utm.to_crs(epsg=4326)
        patch_densities_latlon.to_file(f'{folder_name}geodata/block_densities_latlon.geojson', driver='GeoJSON')
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
        
        blocks_latlon.to_file(folder_name+'geodata/'+'blocks'+'latlon'+'.geojson', driver='GeoJSON')
    
    
    ft = datetime.datetime.now()
    return ft-dt
    

def people_near_x(folder_name, geodata_path, boundaries, year, utm_crs, sqkm_per_pixel):
    if os.path.exists(geodata_path):
        service_latlon = gpd.read_file(geodata_path)
        service_latlon = service_latlon.intersection(boundaries)
        service_mw = service_latlon.to_crs('ESRI:54009')
        service_area = service_latlon.to_crs(utm_crs).area.sum()
        if service_area == 0:
            total_PNS = 0
        else:
            modulo = year % 5
            earlier = year - modulo
            later = year + (5 - modulo)
            earlier_stats = rasterstats.zonal_stats(
                service_mw,
                f"{folder_name}geodata/population/pop_{earlier}.tif", 
                stats=['mean'], 
                all_touched=True
                ) 
            earlier_dens = earlier_stats[0]['mean'] / sqkm_per_pixel 
            later_stats = rasterstats.zonal_stats(
                service_mw,
                f"{folder_name}geodata/population/pop_{later}.tif", 
                stats=['mean'], 
                all_touched=True
                ) 
            later_dens = later_stats[0]['mean'] / sqkm_per_pixel 
            peryear_diff = (later_dens - earlier_dens) / 5
            mean_dens_per_m2 = (earlier_dens + (modulo * peryear_diff)) / 1000000 #km to m
            total_PNS = mean_dens_per_m2 * service_area
    else:
        total_PNS = 0 
    return total_PNS


def calculate_indicators(boundaries, 
                      folder_name='', 
                      to_test = [
                           'healthcare',
                           'schools',
                           'h+s',
                           #'libraries',
                           'bikeshare',
                           'carfree',
                           'blocks',
                           'density',
                           'pnft',
                           'pnrt',
                           'pnpb', #protected bikeways
                           'pnab', #all bikeways
                           'highways',
                           #'special',
                           #'transport_performance',
                           #'connectome'
                           ],
                      #years = range(1975,2031),
                      years = [1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2025],
                      current_year = 2022,
                      ghsl_resolution = '1000',
                      debug = True,
                      ):   

    if folder_name != '' and not folder_name[-1:] == '/':
        folder_name += '/'
        
    results = {}
     
    #TODO change this (and similar above) to osmnx function 
    boundaries_latlon = gpd.GeoDataFrame(geometry=[boundaries], crs=4326)
    longitude = round(numpy.mean(boundaries_latlon.geometry.centroid.x),10)
    utm_zone = int(math.floor((longitude + 180) / 6) + 1)
    utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
    boundaries_utm = boundaries_latlon.to_crs(utm_crs)
    boundaries_mw = boundaries_utm.to_crs("ESRI:54009")
    sqkm_per_pixel = (float(ghsl_resolution) / 1000) ** 2
        
    total_pops = {}
    for year in years:
        if year % 5 == 0:
            pop_stats = rasterstats.zonal_stats(
                boundaries_mw,
                f"{folder_name}geodata/population/pop_{year}.tif", 
                stats=['mean'], 
                all_touched=True
                ) 
            mean_density_per_km2 = pop_stats[0]['mean'] / sqkm_per_pixel
            mean_density_per_m2 = mean_density_per_km2 / 1000000
            total_pop = mean_density_per_m2 * boundaries_utm.area.sum()
            print(f'Total pop {year}:', total_pop)
            total_pops[year] = total_pop
            
    for year in years:
        modulo = year % 5
        if modulo == 0:
            results[f'total_pop_{year}'] = total_pops[year]
        else:
            earlier = year - modulo
            later = year + (5 - modulo)
            earlier_pop = total_pops[earlier] 
            later_pop = total_pops[later] 
            peryear_diff = (later_pop - earlier_pop) / 5
            total_pops[year] = earlier_pop + (modulo * peryear_diff)
            results[f'total_pop_{year}'] = earlier_pop + (modulo * peryear_diff)
            
    if total_pops[current_year] == 0:
        results['recording_time'] = str(datetime.datetime.now())
        return results
        
    if 'density' in to_test:
        densities = {}
        for year in years:
            if year % 5 == 0:
                density = rasterstats.zonal_stats(boundaries_mw, 
                                        f"{folder_name}geodata/population/pop_{year}.tif", 
                                        stats = [],
                                        add_stats={'weighted': weighted_pop_density}
                                        )[0]['weighted']
                densities[year] = density / sqkm_per_pixel 
                print('weighted pop density', year, density / sqkm_per_pixel)
        for year in years:
            modulo = year % 5
            if modulo == 0:
                results[f'density_{year}'] = densities[year]
            else:
                earlier = year - modulo
                later = year + (5 - modulo)
                earlier_dens = densities[earlier]
                later_dens = densities[later] 
                peryear_diff = (later_dens - earlier_dens) / 5
                results[f'density_{year}'] = earlier_dens + (modulo * peryear_diff)
    
    services = ['healthcare','schools','h+s','libraries','bikeshare','pnab','pnpb',
                'pnft','carfree','special']
    for service in services:
        if service in to_test:
            geodata_path = f'{folder_name}geodata/{service}latlon.geojson'
            total_PNS = people_near_x(folder_name, geodata_path, boundaries, current_year, utm_crs, sqkm_per_pixel)
            print('Total People Near Service for', service, ":", total_PNS, 100*total_PNS/total_pops[current_year],"%")
            results[f'{service}_{current_year}'] = total_PNS / total_pops[current_year]
            
    for service_with_points in ['healthcare', 'schools', 'libraries', 'bikeshare', 'pnft','special',]:
        if service_with_points in to_test:
            geodata_path = folder_name+'geodata/'+service_with_points+'_points_latlon'+'.geojson'
            total_services = gpd.read_file(geodata_path).intersection(boundaries)
            results[f'n_points_{service}_{current_year}'] = len(total_services)
        
            
    if 'pnab' in to_test:
        if os.path.exists(folder_name+'geodata/allbike_latlon.geojson'):
            all_bikeways_latlon = gpd.read_file(folder_name+'geodata/allbike_latlon.geojson')
            all_bikeways_latlon = all_bikeways_latlon.intersection(boundaries)
            unprotected_m = sum(all_bikeways_latlon.to_crs(utm_crs).geometry.length)
        else:
            unprotected_m = 0
        results['all_bikeways_km'] = unprotected_m / 1000
        
    if 'pnpb' in to_test:
        if os.path.exists(folder_name+'geodata/protectedbike_latlon.geojson'):
            protected_bikeways_latlon = gpd.read_file(f'{folder_name}geodata/protectedbike_latlon.geojson')
            protected_bikeways_latlon = protected_bikeways_latlon.intersection(boundaries)
            protected_m = sum(protected_bikeways_latlon.to_crs(utm_crs).geometry.length)
        else:
            protected_m = 0
        results['protected_bikeways_km'] = protected_m / 1000
        
    if 'highways' in to_test:
        geodata_path = f'{folder_name}geodata/buffered_hwys_latlon.geojson'
        near_hwys = people_near_x(folder_name, geodata_path, boundaries, current_year, utm_crs, sqkm_per_pixel)
        not_near_hwys = total_pops[current_year] - near_hwys
        print('Total People Safe From Highways:', not_near_hwys, 100*not_near_hwys/total_pops[current_year],"%")
        results['people_not_near_highways'] = not_near_hwys / total_pops[current_year]
        
        if os.path.exists(folder_name+'geodata/allhwys_latlon.geojson'):
            hwypoly_latlon = gpd.read_file(f'{folder_name}geodata/allhwys_latlon.geojson')
            hwypoly_latlon = hwypoly_latlon.intersection(boundaries)
            hwy_m = sum(protected_bikeways_latlon.to_crs(utm_crs).geometry.length) / 4 #divide by 4 because we're looking at divided highway polys, not lines :)
        else:
            hwy_m = 0
        results['highway_km'] = hwy_m / 1000
    
    if 'pnrt' in to_test: 
        geodata_path = f'{folder_name}geodata/rapid_transit/{current_year}/all_isochrones_ll.geojson'
        if os.path.exists(geodata_path):
            rt_gdf = gpd.read_file(geodata_path)
            rt_in_bounds = rt_gdf.intersection(boundaries)
            if rt_in_bounds.unary_union is not None:
                for year in years:
                    if year <= current_year: #TODO efficiency - intersections BEFORE mode loop
                        for mode in ['all','mrt','lrt','brt']:
                            #PNRT
                            geodata_path = f'{folder_name}geodata/rapid_transit/{year}/{mode}_isochrones_ll.geojson'
                            total_pnrt = people_near_x(folder_name, geodata_path, boundaries, year, utm_crs, sqkm_per_pixel)
                            if mode == "all" and year % 5 == 0:
                                print (f"all-modes PNrT {year}: {100*total_pnrt/total_pops[year]}")
                            results[f'PNrT_{mode}_{year}'] = total_pnrt/total_pops[year]
                            # # stations, kms of line, RTR
                            if os.path.exists(f'{folder_name}geodata/rapid_transit/{year}/{mode}_lines_ll.geojson'):
                                lines_ll = gpd.read_file(f'{folder_name}geodata/rapid_transit/{year}/{mode}_lines_ll.geojson')
                                lines_ll = lines_ll.intersection(boundaries)
                                km = sum(lines_ll.to_crs(utm_crs).geometry.length) / 1000
                            else:
                                km=0
                            if os.path.exists(f'{folder_name}geodata/rapid_transit/{year}/{mode}_stations_ll.geojson'):
                                stns_ll = gpd.read_file(f'{folder_name}geodata/rapid_transit/{year}/{mode}_stations_ll.geojson')
                                stns_ll = stns_ll.intersection(boundaries)
                                n_stns = len(stns_ll)
                            else:
                                n_stns = 0
                            results[f'km_{mode}_{year}'] = km
                            results[f'stns_{mode}_{year}'] = n_stns
                            results[f'rtr_{mode}_{year}'] = km / (total_pops[year]/1000000)
            else:
                for year in years:
                    if year <= current_year:
                        for mode in ['all','mrt','lrt','brt']:
                            results[f'PNrT_{mode}_{year}'] = 0
                            results[f'km_{mode}_{year}'] = 0
                            results[f'stns_{mode}_{year}'] = 0
                            results[f'rtr_{mode}_{year}'] = 0
        else:
            for year in years:
                if year <= current_year:
                    for mode in ['all','mrt','lrt','brt']:
                        results[f'PNrT_{mode}_{year}'] = 0
                        results[f'km_{mode}_{year}'] = 0
                        results[f'stns_{mode}_{year}'] = 0
                        results[f'rtr_{mode}_{year}'] = 0
                        
    if 'blocks' in to_test:
        geodata_path = folder_name+'geodata/'+'blocks'+'latlon'+'.geojson'
        if os.path.exists(geodata_path):
            blocks = gpd.read_file(geodata_path)
            selection = blocks[blocks.intersects(boundaries)]
            av_size = selection.area_utm.mean()
            block_density = 1000000 / av_size
        else:
            block_density = 'NA'
        results['block_density'] = block_density
            
            
            
    if 'transport_performance' in to_test:
        # any_gtfs = False
        # for file in os.listdir(folder_name+'temp/access/'):
        #     if file[-4] == '.zip':
        #         any_gtfs = True
        cxn_gdf_latlon = gpd.read_file(f'{folder_name}geodata/connections.gpkg')
        cxn_gdf_overlap = cxn_gdf_latlon.overlay(boundaries_latlon, how='intersection')
        for column in cxn_gdf_overlap.columns:
            name_parts = column.split('_')
            if name_parts[0] == 'performance':
                total_val = sum(cxn_gdf_overlap[column] * cxn_gdf_overlap['population'])
                weighted_avg = total_val/sum(cxn_gdf_overlap['population'])
                print(column, ': ',weighted_avg)
                results[column] = weighted_avg
        
            
    results['recording_time'] = str(datetime.datetime.now())
        
    return results

