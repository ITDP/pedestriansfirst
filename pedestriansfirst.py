import fiona
import datetime
import os
import numpy
import math
import statistics
import rasterstats
import rasterio.mask
import subprocess
import json
import traceback
import utm
import geojson
import shutil
#import utm_zone

import numpy as np
import osmnx as ox
import networkx as nx
import geopandas as gpd
import shapely.geometry
from shapely.geometry import LineString
from shapely.ops import unary_union
import shapely.ops

import isochrones
import get_service_locations
import gtfs_parser

import pdb
import logging


ox.utils.config(log_console = False)

def make_patches(boundaries, patch_length = 5): #patch_length in km
    ''' 
    'tile' the boundaries of a city into patches, like a patchwork quilt
    '''
    
    longitude_factor = 0.00898 # degrees per km
    longitude_factor_m = 0.00898 / 1000 # degrees per m
    latitude_factor = (math.cos(abs(boundaries.bounds[1])*0.0174533))/111.319
    latitude_factor_m = latitude_factor / 1000
    
    bbox = boundaries.bounds
    
    height_degrees = abs(bbox[3]-bbox[1])
    height_km = height_degrees / latitude_factor
    width_degrees = abs(bbox[2]-bbox[0])
    width_km = width_degrees / longitude_factor
    
    n_hslicers = math.floor(height_km / patch_length)
    n_vslicers = math.floor(width_km / patch_length)
    
    hslicers = []
    vslicers = []
    
    for i in range(1,n_hslicers+1):
        h_increment = (bbox[3]-bbox[1])/(n_hslicers+1)
        lat = bbox[1]+(i*h_increment)
        slicer = shapely.geometry.LineString([(bbox[0],lat),(bbox[2],lat)])
        hslicers.append(slicer)
        
    for i in range(1,n_vslicers+1):
        v_increment = (bbox[2]-bbox[0])/(n_vslicers+1)
        lon = bbox[0]+(i*v_increment)
        slicer = shapely.geometry.LineString([(lon,bbox[1]),(lon, bbox[3])])
        hslicers.append(slicer)
        
    patches = shapely.geometry.MultiPolygon(polygons=[boundaries])
    for slicer in hslicers+vslicers:
        patches = shapely.geometry.MultiPolygon(polygons = shapely.ops.split(patches, slicer))
    
    print("cut"+str(len(list(patches)))+"patches")
    
    return list(patches)

def weighted_pop_density(array):
    total = 0
    for cell in array.compressed():
        total += cell**2
    return total / numpy.sum(array)

def pedestrians_first(boundaries, 
                      id_code,
                      name,
                      folder_name='', 
                      buffer_dist=100,#m
                      headway_threshold=10,#min
                      to_test = [
                           'healthcare',
                           'schools',
                           'h+s',
                           'libraries',
                           'bikeshare',
                           'carfree',
                           'blocks',
                           'density',
                           'transit',
                           'pnpb', #protected bikeways
                           'pnab', #all bikeways
                           #'special',
                           ],
                      distances = { #network buffers, in meters
                            'healthcare': 1000,
                            'schools': 1000,
                            'libraries': 1000,
                            'bikeshare': 500,
                            'transit': 500,
                            'special': 500,
                            'pnpb': 250,
                            'pnab': 250,
                            },
                      overpass = False,
                      patch_length = 5, #km
                      boundary_buffer = 0, #km
                      blocks_simplification = 15, #m
                      gtfs_files = [],
                      ):    
    dt = datetime.datetime.now()
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    fh = logging.FileHandler('log_filename.txt')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    useful_tags = ox.settings.useful_tags_way + ['cycleway', 'cycleway:left', 'cycleway:right']
    ox.config(use_cache=True, log_console=True, useful_tags_way=useful_tags)
    
    if folder_name != '' and not folder_name[-1:] == '/':
        folder_name += '/'
    
    bound_latlon = gpd.GeoDataFrame(geometry = [boundaries])
    bound_latlon.crs = {'init':'epsg:4326'}
    longitude = round(numpy.mean(bound_latlon.geometry.centroid.x),10)
    utm_zone = int(math.floor((longitude + 180) / 6) + 1)
    utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
    crs=utm_crs
    
    if boundary_buffer > 0:
        bound_utm = bound_latlon.to_crs(utm_crs)
        bound_utm.geometry = bound_utm.geometry.buffer(boundary_buffer*1000)
        bound_latlon = bound_utm.to_crs(epsg=4326)
        boundaries = bound_latlon.geometry.unary_union
    print(utm_crs)
    bbox = boundaries.bounds
    
    reppoint = boundaries.representative_point()
    #crs = utm_zone.epsg(geojson.Point((reppoint.x[0],reppoint.y[0])))

    
    if type(boundaries) == shapely.geometry.multipolygon.MultiPolygon:
        patches = []
        for poly in list(boundaries):
            patches += make_patches(poly, patch_length=patch_length)
    else:
        patches = make_patches(boundaries, patch_length=patch_length)
    
    longitude_factor = 0.00898 # degrees per km
    longitude_factor_m = 0.00898 / 1000 # degrees per m
    latitude_factor = (math.cos(abs(boundaries.bounds[1])*0.0174533))/111.319
    latitude_factor_m = latitude_factor / 1000
    
    print('Evaluating Pedestrians First indicators in',name)
    print('Measuring',str(to_test))
    
    quilt_isochrone_polys = {}
    quilt_center_nodes = {}
    for service in to_test:
        quilt_isochrone_polys[service] = False
        quilt_center_nodes[service] = []
        if 'pnab' in to_test or 'pnpb' in to_test:
            quilt_isochrone_polys['pnab'] = False
            quilt_center_nodes['pnab'] = []
            quilt_isochrone_polys['pnpb'] = False
            quilt_center_nodes['pnpb'] = []
            
    
    patch_times = []
    
    all_coords={}
    
    testing_services = []
    for service in ['healthcare', 'schools', 'libraries', 'bikeshare']:
        if service in to_test:
            testing_services.append(service)
            #all_coords[service] = get_point_locations(boundaries, queries[service])
    
    if len(testing_services) > 0:
        if overpass:
            raise RuntimeError
            #haven't set this up yet, because I don't have a query for bikeshare or carfree
        else:
            handler = get_service_locations.ServiceHandler()
            handler.apply_file(folder_name+'city.o5m', locations=True)
            for service in testing_services:
                all_coords[service] = handler.locationlist[service]
                citywide_carfree = handler.carfreelist
            
    #should I finish this?
    if 'highways' in to_test:     
        highway_filter = '["area"!~"yes"]["highway"~"motorway|trunk"]'
        G_hwys = ox.graph_from_polygon(boundaries,
                                custom_filter=highway_filter,
                                retain_all=True)
        G_hwys = ox.project_graph(G_hwys)
        edges_hwys = ox.utils_graph.graph_to_gdfs(G_hwys, nodes=False)
        edges_polyline = edges_hwys.geometry.unary_union
        increment_dist = 25
        distances = np.arange(0, edges_polyline.length, increment_dist)
        points = [edges_polyline.interpolate(distance) for distance in distances] + [edges_polyline.boundary[1]]
        multipoint = unary_union(points)
        #need to switch back to latlon and put in all_coords['highway']
        
    if 'special' in to_test:
        if os.path.isfile(folder_name+'special.shp'):
            testing_services.append('special')
            special = gpd.read_file(folder_name+'special.shp')
            all_coords['special'] = [(pt.y, pt.x) for pt in special.geometry]

    if 'transit' in to_test:
        testing_services.append('transit')
        #this approach does not let us combine OpenMobilityData files with local files.
        if len(gtfs_files) == 0:
            sources = gtfs_parser.get_feed_infos(gtfs_parser.get_relevant_locs(bbox))
            #note! This doesn't reflect a buffered boundary if boundary_buffer > 0
            transit_stop_sets = gtfs_parser.count_all_sources(sources, 
                                                              source_type = 'openmobilitydata',
                                                              headwaylim = headway_threshold * 2)
        else:
            transit_stop_sets = gtfs_parser.count_all_sources(gtfs_files, 
                                                              source_type = 'local_files',
                                                              headwaylim = headway_threshold * 2)
            
    if 'pnab' in to_test or 'pnpb' in to_test:
        testing_services.append('pnpb')
        testing_services.append('pnab')
        total_protectedbike_oneway = gpd.GeoDataFrame()
        total_protectedbike_twoway = gpd.GeoDataFrame()
        total_unprotectedbike_oneway = gpd.GeoDataFrame()
        total_unprotectedbike_twoway = gpd.GeoDataFrame()
        
    if len(to_test) > 0 and to_test != ["blocks"]:
        for p_idx, patch in enumerate(patches):
            try:
                
            
                patch_start = datetime.datetime.now()
                
                unbuffered_patch = patch
                
                max_service_dist_km = max(distances.values())/1000
                
                # patch = shapely.geometry.box(
                        # patch.bounds[0] - (max_service_dist_km * longitude_factor),
                        # patch.bounds[1] - (max_service_dist_km * latitude_factor),
                        # patch.bounds[2] + (max_service_dist_km * longitude_factor),
                        # patch.bounds[3] + (max_service_dist_km * latitude_factor)
                        # )
                        
                patchgdf = gpd.GeoDataFrame(geometry=[patch], crs=4326)
                patchgdf_utm = patchgdf.to_crs(crs) #NEED TO DEFINE CRS EARLIER
                patchgdf_utm.geometry = patchgdf_utm.geometry.buffer(max_service_dist_km*1000)
                patchgdf = patchgdf_utm.to_crs(4326)
                patch = patchgdf.geometry.unary_union.intersection(boundaries)
                patchgdf = gpd.GeoDataFrame(geometry=[patch], crs=4326)
                
                if os.path.exists('pathbounds.geojson'):
                    os.remove('patchbounds.geojson')
                if os.path.exists('pathbounds.poly'):
                    os.remove('patchbounds.poly')
                if os.path.exists('patch.osm'):
                    os.remove('patch.osm')
                if os.path.exists('allhwyspatch.osm'):
                    os.remove('allhwyspatch.osm')
                
                patchgdf.to_file('patchbounds.geojson', driver='GeoJSON')
                
                if os.path.exists('patchbounds.poly'):
                    os.remove('patchbounds.poly')
                shutil.copy('patchbounds.geojson', f'{str(folder_name)}/patchbounds{str(p_idx)}.geojson')
                subprocess.run('python ogr2poly/ogr2poly.py patchbounds.geojson > patchbounds.poly', shell=True, check=True)
                
                walk_filter = ('["area"!~"yes"]["highway"!~"link|motor'
                               '|proposed|construction|abandoned'
                               '|platform|raceway"]'
                               '["service"!~"parking_aisle|driveway"]'
                               '["foot"!~"no"]["service"!~"private"]'
                               '{}').format(ox.settings.default_access)
                
                if overpass:
                    G = ox.graph_from_polygon(patch, 
                                              custom_filter=walk_filter, 
                                              simplify=False, 
                                              retain_all=True)
                    G_allhwys = G.copy()
                else:
                    try:
                        boundingarg = '-b='
                        boundingarg += str(patch.bounds[0])+','
                        boundingarg += str(patch.bounds[1])+','
                        boundingarg += str(patch.bounds[2])+','
                        boundingarg += str(patch.bounds[3])
                        subprocess.check_call(['osmconvert',
                                               str(folder_name)+'citywalk.o5m',
                                               #boundingarg, #OR
                                               "-B=patchbounds.poly",
                                               #'--complete-ways',  #was commented
                                               '--drop-broken-refs',  #was uncommented
                                               '-o=patch.osm'])
                        G = ox.graph_from_xml('patch.osm', 
                                               simplify=False, retain_all=True)
                        os.remove('patch.osm')
                        if 'pnab' in to_test or 'pnpb' in to_test:
                            subprocess.check_call(['osmconvert',
                                               str(folder_name)+'cityhighways.o5m',
                                               #boundingarg, #OR
                                               "-B=patchbounds.poly",
                                               #'--complete-ways',
                                               '--drop-broken-refs',
                                               '-o=allhwyspatch.osm'])
                            G_allhwys = ox.graph_from_xml('allhwyspatch.osm', 
                                               simplify=False, retain_all=True)
                            os.remove('allhwyspatch.osm')
                    except KeyError: #something to do with clipping, seems to happen once in a while
                        #this is a very stupid band-aid, but it works for now, I think
                        print ('KEYERROR FROM CLIPPING PATCH', p_idx)
                        with open(str(folder_name)+"patcherrorlog.txt", "a") as patcherrorlog:
                            patcherrorlog.write('KEYERROR FROM CLIPPING PATCH '+str(p_idx))
                        G = ox.graph_from_polygon(patch, 
                                              custom_filter=walk_filter, 
                                              simplify=False, 
                                              retain_all=True)
                        G_allhwys = G.copy()
                os.remove('patchbounds.poly')
                        
                
                G.remove_nodes_from(list(nx.isolates(G)))
                if 'pnab' in to_test or 'pnpb' in to_test:
                    G_allhwys.remove_nodes_from(list(nx.isolates(G_allhwys)))
                
                simple_G = ox.simplify_graph(G)
                
                center_nodes = {}
                for service in all_coords.keys():
                    if service in ['healthcare','schools','libraries','bikeshare','special']:
                        center_nodes[service] = []
                        for coord in all_coords[service]:
                            lat = coord[0]
                            lon = coord[1]
                            if unbuffered_patch.bounds[0] < lon < unbuffered_patch.bounds[2] and unbuffered_patch.bounds[1] < lat < unbuffered_patch.bounds[3]:
                                point = shapely.geometry.Point(lon,lat)
                                nearest = ox.get_nearest_node(simple_G, coord)
                                if not nearest in center_nodes[service]:    
                                    center_nodes[service].append(nearest)
                
                if 'pnab' in to_test or 'pnpb' in to_test:
                    allhwys_gdf = ox.graph_to_gdfs(G_allhwys, nodes=False)
                    print(allhwys_gdf.columns)
                    protected_oneway = gpd.GeoDataFrame()
                    protected_twoway = gpd.GeoDataFrame()
                    unprotected_oneway = gpd.GeoDataFrame()
                    unprotected_twoway = gpd.GeoDataFrame()
                    if 'highway' in allhwys_gdf.columns:
                        if 'oneway' in allhwys_gdf.columns:
                            protected_oneway = protected_oneway.append(allhwys_gdf[(allhwys_gdf['highway'] == 'cycleway') & (allhwys_gdf['oneway'] == 'yes')])
                            protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['highway'] == 'cycleway') & (allhwys_gdf['oneway'] != 'yes')])
                        else:
                            protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['highway'] == 'cycleway')])
                    if ('highway' in allhwys_gdf.columns) and ('bicycle' in allhwys_gdf.columns):
                        if 'oneway' in allhwys_gdf.columns:
                            protected_oneway = protected_oneway.append(allhwys_gdf[(allhwys_gdf['highway'] == 'path') & (allhwys_gdf['bicycle'] == 'designated') & (allhwys_gdf['oneway'] == 'yes')])
                            protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['highway'] == 'path') & (allhwys_gdf['bicycle'] == 'designated') & (allhwys_gdf['oneway'] != 'yes')])
                        else:
                            protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['highway'] == 'path') & (allhwys_gdf['bicycle'] == 'designated')])
                    if 'cycleway' in allhwys_gdf.columns:
                        if 'oneway' in allhwys_gdf.columns:
                            if 'oneway:bicycle' in allhwys_gdf.columns:
                                protected_oneway = protected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'track') & (allhwys_gdf['oneway'] == 'yes') & (allhwys_gdf['oneway:bicycle'] != 'no')])
                                protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'track') & (allhwys_gdf['oneway'] == 'yes') & (allhwys_gdf['oneway:bicycle'] == 'no')])
                                protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'track') & (allhwys_gdf['oneway'] != 'yes')])
                                unprotected_oneway = unprotected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'lane') & (allhwys_gdf['oneway'] == 'yes') & (allhwys_gdf['oneway:bicycle'] != 'no')])
                                unprotected_twoway = unprotected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'lane') & (allhwys_gdf['oneway'] == 'yes') & (allhwys_gdf['oneway:bicycle'] == 'no')])
                                unprotected_twoway = unprotected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'lane') & (allhwys_gdf['oneway'] != 'yes')])
                            else:
                                protected_oneway = protected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'track') & (allhwys_gdf['oneway'] == 'yes')])
                                protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'track') & (allhwys_gdf['oneway'] != 'yes')])
                                unprotected_oneway = unprotected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'lane') & (allhwys_gdf['oneway'] == 'yes')])
                                unprotected_twoway = unprotected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'lane') & (allhwys_gdf['oneway'] != 'yes')])
                        else:
                            protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'track')])
                            unprotected_twoway = unprotected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway'] == 'lane')])
                    if 'cycleway:both' in allhwys_gdf.columns:
                        protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway:both'] == 'track')])
                        unprotected_twoway = unprotected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway:both'] == 'lane')])
                    if 'cycleway:left' in allhwys_gdf.columns:
                        if 'cycleway:left:oneway' in allhwys_gdf.columns:
                            protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway:left'] == 'track') & (allhwys_gdf['cycleway:left:oneway'] == 'no')])
                            unprotected_twoway = unprotected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:left'] == 'lane') & (allhwys_gdf['cycleway:left:oneway'] == 'no')])
                            protected_oneway = protected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:left'] == 'track') & (allhwys_gdf['cycleway:left:oneway'] != 'no')])
                            unprotected_oneway = unprotected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:left'] == 'lane') & (allhwys_gdf['cycleway:left:oneway'] != 'no')])
                        else:
                            protected_oneway = protected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:left'] == 'track')])
                            unprotected_oneway = unprotected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:left'] == 'lane')])
                    if 'cycleway:right' in allhwys_gdf.columns:
                        if 'cycleway:right:oneway' in allhwys_gdf.columns:
                            protected_twoway = protected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway:right'] == 'track') & (allhwys_gdf['cycleway:right:oneway'] == 'no')])
                            unprotected_twoway = unprotected_twoway.append(allhwys_gdf[(allhwys_gdf['cycleway:right'] == 'lane') & (allhwys_gdf['cycleway:right:oneway'] == 'no')])
                            protected_oneway = protected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:right'] == 'track') & (allhwys_gdf['cycleway:right:oneway'] != 'no')])
                            unprotected_oneway = unprotected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:right'] == 'lane') & (allhwys_gdf['cycleway:right:oneway'] != 'no')])
                        else:
                            protected_oneway = protected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:right'] == 'track')])
                            unprotected_oneway = unprotected_oneway.append(allhwys_gdf[(allhwys_gdf['cycleway:right'] == 'lane')])
                    total_protectedbike_oneway = total_protectedbike_oneway.append(protected_oneway)
                    total_unprotectedbike_oneway = total_unprotectedbike_oneway.append(unprotected_oneway)
                    total_protectedbike_twoway = total_protectedbike_twoway.append(protected_twoway)
                    total_unprotectedbike_twoway = total_unprotectedbike_twoway.append(unprotected_twoway)
                    center_nodes['pnpb'] = []
                    center_nodes['pnab'] = []
                    pnpbnodes = []
                    pnabnodes = []
                    for edge in protected_twoway.index:
                        if not edge[0] in pnpbnodes:
                            pnpbnodes.append(edge[0])
                        if not edge[1] in pnpbnodes:
                            pnpbnodes.append(edge[1])
                        if not edge[0] in pnabnodes:
                            pnabnodes.append(edge[0])
                        if not edge[1] in pnabnodes:
                            pnabnodes.append(edge[1])
                    for edge in protected_oneway.index:
                        if not edge[0] in pnpbnodes:
                            pnpbnodes.append(edge[0])
                        if not edge[1] in pnpbnodes:
                            pnpbnodes.append(edge[1])
                        if not edge[0] in pnabnodes:
                            pnabnodes.append(edge[0])
                        if not edge[1] in pnabnodes:
                            pnabnodes.append(edge[1])
                    for edge in unprotected_oneway.index:
                        if not edge[0] in pnabnodes:
                            pnabnodes.append(edge[0])
                        if not edge[1] in pnabnodes:
                            pnabnodes.append(edge[1])
                    for edge in unprotected_twoway.index:
                        if not edge[0] in pnabnodes:
                            pnabnodes.append(edge[0])
                        if not edge[1] in pnabnodes:
                            pnabnodes.append(edge[1])
                    for node in pnpbnodes:
                        if node in simple_G.nodes:
                            center_nodes['pnpb'].append(node)
                    for node in pnabnodes:
                        if node in simple_G.nodes:
                            center_nodes['pnab'].append(node)
                
                if 'transit' in to_test:
                    center_nodes['transit'] = []
                    transit_centers = {}
                    for service_idx, service in enumerate(transit_stop_sets):
                        for stop_id in service.keys():
                            lat = float(service[stop_id][1])
                            lon = float(service[stop_id][2])
                            headway = float(service[stop_id][0])
                            if unbuffered_patch.bounds[0] < lon < unbuffered_patch.bounds[2] and unbuffered_patch.bounds[1] < lat < unbuffered_patch.bounds[3]:
                                center_node = ox.get_nearest_node(simple_G, (lat, lon))
                                if center_node not in transit_centers:
                                    transit_centers[center_node] = {service_idx : headway} #store the headway value
                                elif service_idx not in transit_centers[center_node].keys():
                                    transit_centers[center_node][service_idx] = headway
                                elif headway < transit_centers[center_node][service_idx]:
                                    transit_centers[center_node][service_idx] = headway
                    for center_node in transit_centers.keys():
                        if len(transit_centers[center_node]) > 0:
                            inv_headway = 0
                            min_hw = 100
                            for service_idx in transit_centers[center_node].keys():
                                inv_headway += 1 / transit_centers[center_node][service_idx]
                                if transit_centers[center_node][service_idx] < min_hw:
                                    min_hw = transit_centers[center_node][service_idx]
                            headway = 1 / inv_headway
                            if headway <= headway_threshold:
                                center_nodes['transit'].append(center_node)
                                
                
                
                
                
                # Project Graph
                if not crs: #this only runs once, crs is invariable over patches
                    G = ox.project_graph(G)
                    crs = G.graph['crs'] 
                else:
                    G = ox.project_graph(G, to_crs=crs)
                G = ox.simplify_graph(G)
                
                isochrone_polys = {}
                failures = {}
                for service in to_test:
                    failures[service] = 0
                    if 'pnab' in to_test or 'pnpb' in to_test:
                        failures['pnab']=0
                        failures['pnpb']=0
                
                
                
                # Get polygons
                for service in testing_services:
                    print('getting polygons for',service,len(center_nodes[service]),'center_nodes')
                    isochrone_polys[service], fails = isochrones.make_iso_polys(G, center_nodes[service], distance=distances[service], edge_buff=buffer_dist)
                    failures[service] += fails
                    
                for service in isochrone_polys.keys():
                    if service not in quilt_isochrone_polys.keys() or not quilt_isochrone_polys[service]:
                        quilt_isochrone_polys[service] = isochrone_polys[service]
                    elif isochrone_polys[service]:
                        quilt_isochrone_polys[service] = shapely.ops.cascaded_union([quilt_isochrone_polys[service],isochrone_polys[service]])
                for service in center_nodes.keys():
                    quilt_center_nodes[service] = quilt_center_nodes[service] + center_nodes[service]
                
                patch_time = datetime.datetime.now() - patch_start
                    
                print("finished patch #",p_idx,'out of', len(patches),"in",str(patch_time))
                patch_times.append(patch_time)
            except ox._errors.EmptyOverpassResponse:
                print('EmptyOverpassResponse')
                pass
            except ValueError:
                print('ValueError')
                pass #sorry
            #except:
            #    print("GOT SOME ERROR FOR PATCH", p_idx)
            # except ValueError:
            #     print('ValueError')
            #     now = str(datetime.datetime.now())
            #     with open('error'+now+'.txt','w') as errout:
            #         traceback.print_exc(limit=3,file=errout)
            #     print('saved to error'+now+'.txt')
    #epsg = 32600+int(crs.split(' ')[1].split('=')[1]) #This is wild -- 
    #osmnx seems to just give all data in northern-hemisphere format
    #Sorry about the stupid parsing of the projection definition, I'm lazy 
    
    results = {'name':name,'id_code':id_code}
    
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
     
    boundaries_latlon = gpd.GeoDataFrame(geometry=[boundaries])
    boundaries_latlon.crs = {'init':'epsg:4326'}
    try:
        boundaries_utm = boundaries_latlon.to_crs(crs)
    except ValueError:
        longitude = round(numpy.mean(boundaries_latlon.geometry.centroid.x),10)
        utm_zone = int(math.floor((longitude + 180) / 6) + 1)
        utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
        crs = utm_crs
        boundaries_utm = boundaries_latlon.to_crs(crs)
        
    stats = rasterstats.zonal_stats(boundaries_latlon, 'pop_dens.tif', stats=['sum'])           
    total_pop = stats[0]['sum'] 
    
    for service in testing_services:
        if quilt_isochrone_polys[service]:
            service_utm = gpd.GeoDataFrame(geometry = [quilt_isochrone_polys[service]])
            service_utm.crs = crs#{'init':'epsg:'+str(epsg)}
            service_utm.geometry = service_utm.geometry.simplify(15) #maybe this should be after the population calculation
            service_utm = gpd.overlay(service_utm ,boundaries_utm, how='intersection')
            service_latlon = service_utm.to_crs(epsg=4326)
            service_latlon.to_file(folder_name+service+'latlon'+'.geojson', driver='GeoJSON')
            
            stats = rasterstats.zonal_stats(service_latlon, 'pop_dens.tif', stats=['sum'])           
            total_PNS = stats[0]['sum']
            print("\n")
            print('Total People Near Service for', service, ":", total_PNS)
            print(100*total_PNS/total_pop,"% of",total_pop)
            results[service] = total_PNS / total_pop
        else:
            print ('NO SERVICE FOR', service)
            results[service] = 0
            
    if 'pnpb' in to_test or 'pnab' in to_test:
        protected_km = 0
        unprotected_km = 0
        if not os.path.exists(folder_name+'bikeways/'):
            os.makedirs(folder_name+'bikeways/')
        if not total_protectedbike_oneway.empty:
            merged_protectedbike_oneway = gpd.GeoDataFrame(geometry = [total_protectedbike_oneway.geometry.unary_union], crs=4326)
            merged_protectedbike_oneway.index=range(0,len(merged_protectedbike_oneway))
            merged_protectedbike_oneway = gpd.clip(merged_protectedbike_oneway, boundaries)
            merged_protectedbike_oneway.to_file(folder_name+'bikeways/protected_oneway.geojson',driver='GeoJSON')
            protected_km += sum(merged_protectedbike_oneway.to_crs(crs).geometry.length) / 2 
        if not total_unprotectedbike_oneway.empty:
            merged_unprotectedbike_oneway = gpd.GeoDataFrame(geometry = [total_unprotectedbike_oneway.geometry.unary_union], crs=4326)
            merged_unprotectedbike_oneway = gpd.clip(merged_unprotectedbike_oneway, boundaries)
            merged_unprotectedbike_oneway.crs = 4326
            merged_unprotectedbike_oneway.to_file(folder_name+'bikeways/unprotected_oneway.geojson',driver='GeoJSON')
            unprotected_km += sum(merged_unprotectedbike_oneway.to_crs(crs).geometry.length) / 2 
        if not total_protectedbike_twoway.empty:
            merged_protectedbike_twoway = gpd.GeoDataFrame(geometry = [total_protectedbike_twoway.geometry.unary_union], crs=4326)
            merged_protectedbike_twoway = gpd.clip(merged_protectedbike_twoway, boundaries)
            merged_protectedbike_twoway.crs = 4326
            merged_protectedbike_twoway.to_file(folder_name+'bikeways/protected_twoway.geojson',driver='GeoJSON')
            protected_km += sum(merged_protectedbike_twoway.to_crs(crs).geometry.length)
        if not total_unprotectedbike_twoway.empty:
            merged_unprotectedbike_twoway = gpd.GeoDataFrame(geometry = [total_unprotectedbike_twoway.geometry.unary_union], crs=4326)
            merged_unprotectedbike_twoway = gpd.clip(merged_unprotectedbike_twoway, boundaries)
            merged_unprotectedbike_twoway.crs = 4326
            merged_unprotectedbike_twoway.to_file(folder_name+'bikeways/unprotected_twoway.geojson',driver='GeoJSON')
            protected_km += sum(merged_unprotectedbike_twoway.to_crs(crs).geometry.length) 
        results['protected_bikeways_km'] = protected_km
        results['all_bikeways_km'] = protected_km + unprotected_km
    
    if 'carfree' in to_test:
        print("getting carfree")
        if citywide_carfree:
            carfree_latlon = gpd.GeoDataFrame(geometry = citywide_carfree)
            #just a latlon list of points
            carfree_latlon.crs = {'init':'epsg:4326'}
            carfree_utm = carfree_latlon.to_crs(crs)
            carfree_utm.geometry = carfree_utm.geometry.buffer(100)
            #this is the analysis, the 100m buffer
            carfree_utm = gpd.GeoDataFrame(geometry = [shapely.ops.cascaded_union(carfree_utm.geometry)])
            carfree_utm.geometry = carfree_utm.geometry.simplify(10)
            carfree_utm = gpd.overlay(carfree_utm ,boundaries_utm, how='intersection')
            carfree_utm.crs = crs
            carfree_latlon = carfree_utm.to_crs('epsg:4326')
            
            stats = rasterstats.zonal_stats(carfree_latlon, 'pop_dens.tif', stats=['sum'])
            total_carfree = stats[0]['sum']
            print("\n")
            print('Total People Near Service for carfree', ":", total_carfree)
            print(100*total_carfree/total_pop,"% of",total_pop)
            results['carfree'] = total_carfree / total_pop
            
            carfree_latlon = carfree_utm.to_crs('epsg:4326')
            carfree_latlon.to_file(folder_name+'carfreelatlon'+'.geojson', driver='GeoJSON')
        else:
            print ('NO SERVICE FOR carfree')
            results['carfree'] = 0
            
    
    if 'h+s' in to_test:
        if quilt_isochrone_polys['healthcare'] and quilt_isochrone_polys['schools']:
            service = 'h+s'
            intersect = quilt_isochrone_polys['healthcare'].intersection(quilt_isochrone_polys['schools'])
            if type(intersect) == shapely.geometry.collection.GeometryCollection:
                intersect = [obj for obj in intersect if type(obj) == shapely.geometry.polygon.Polygon]
                intersect = shapely.geometry.MultiPolygon(intersect)
            hs_utm = gpd.GeoDataFrame(geometry = [intersect])
            hs_utm.crs = crs#{'init':'epsg:'+str(epsg)} #maybe this should be after the population calculation
            if hs_utm.geometry.area.sum() != 0:
                hs_utm = gpd.overlay(hs_utm ,boundaries_utm, how='intersection')
                hs_utm.geometry = hs_utm.geometry.simplify(15)
                hs_latlon = hs_utm.to_crs(epsg=4326)
                hs_latlon.to_file(folder_name+service+'latlon'+'.geojson', driver='GeoJSON')
                stats = rasterstats.zonal_stats(hs_latlon, 'pop_dens.tif', stats=['sum'])
                
                total_PNS = stats[0]['sum']
                print("\n")
                print('Total People Near Service for', service, ":", total_PNS)
                print(100*total_PNS/total_pop,"% of",total_pop)
                results[service] = total_PNS / total_pop
            else:
                print ('NO SERVICE FOR h+s')
                results['h+s'] = 0
        else:
            print ('NO SERVICE FOR h+s')
            results['h+s'] = 0
    
    if 'density' in to_test:
        density = rasterstats.zonal_stats(boundaries, 
                                'pop_dens.tif', 
                                stats = [],
                                add_stats={'weighted': weighted_pop_density}
                                )[0]['weighted']
        results['density'] = density / 0.0625 #km^2 / pixel
        print('weighted pop density', results['density'])
        with rasterio.open('pop_dens.tif') as dataset:
            out_image, out_transform = rasterio.mask.mask(dataset, [boundaries], crop=True)
            out_meta = dataset.meta
            out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})
        with rasterio.open(folder_name+"pop_dens.tif", "w", **out_meta) as dest:
            dest.write(out_image)
            
    #garbage collection
    quilt_isochrone_polys = None
    isochrone_polys = None
    out_image = None
    dest = None
    dataset = None
    a = None
    b = None
    import gc 
    gc.collect()
    
    if 'blocks' in to_test:
        print("getting blocks")
        
        patches = make_patches(boundaries, patch_length=patch_length)
        print ("cut", len(list(patches)),"patches for block size in",name)
        
        outblocks = []
        block_counts = []
        for n, patch in enumerate(patches):
            print("patch"+str(n)+" of "+str(len(patches)) )
            unbuffered_patch = gpd.GeoSeries(patch, crs={'init':'epsg:4326'})
            unbuffered_patch = unbuffered_patch.to_crs(crs)[0]
            
            
            buffer_dist = 1
            
            patch = shapely.geometry.box(
                    patch.bounds[0] - (buffer_dist * longitude_factor),
                    patch.bounds[1] - (buffer_dist * latitude_factor),
                    patch.bounds[2] + (buffer_dist * longitude_factor),
                    patch.bounds[3] + (buffer_dist * latitude_factor)
                    )
            
            if overpass:
                G = ox.graph_from_polygon(patch, custom_filter=walk_filter, simplify=False, retain_all=True)
            else:
                boundingarg = '-b='
                boundingarg += str(patch.bounds[0])+','
                boundingarg += str(patch.bounds[1])+','
                boundingarg += str(patch.bounds[2])+','
                boundingarg += str(patch.bounds[3])
                subprocess.check_call(['osmconvert',
                                       folder_name+'citywalk.o5m',
                                       boundingarg,
                                       #'--complete-ways',
                                       '--drop-broken-refs',
                                       '-o=patch.osm'])
                try:
                    G = ox.graph_from_xml('patch.osm', simplify=False, retain_all=True)
                except:
                    G = False
            
            if G:
                try:
                    G = ox.project_graph(G, to_crs=crs)
                    G = ox.simplify_graph(G)
                    
                    streets = ox.utils_graph.graph_to_gdfs(G, nodes = False)
                    
                    if not streets.empty:
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
                                    if block.centroid.within(unbuffered_patch):
                                        area = round(block.area, 3)
                                        perim = round(block.length, 3)
                                        lemgth = round((perim * perim) / area, 3)
                                        if blocks_simplification:
                                            block = block.simplify(blocks_simplification)
                                        all_blocks.append((block, area, perim, lemgth))
                                        if (lemgth < 50) and (1000 < area < 1000000):
                                            selected_areas.append(area)
                            outblocks += all_blocks
                            block_counts.append(len(all_blocks))
                        else:
                            block_counts.append(0)
                            print('not merged!')
                    else:
                        block_counts.append(0)
                except:
                    print('Hawassa Error')
                    block_counts.append(0)
            else:
                block_counts.append(0)
        
        #export            
        
        patch_densities = gpd.GeoDataFrame(geometry = list(patches))
        patch_densities['block_count'] = block_counts
        patch_densities.crs = {'init':'epsg:4326'}
        patch_densities_utm = patch_densities.to_crs(crs)
        patch_densities_utm['density'] = patch_densities_utm.block_count / (patch_densities_utm.area /1000000)
        patch_densities_latlon = patch_densities_utm.to_crs(epsg=4326)
        patch_densities_latlon.to_file(folder_name+'patch_densities'+'latlon'+'.geojson', driver='GeoJSON')
        
        a = gpd.GeoDataFrame(geometry=[block[0] for block in outblocks])
        a.crs = crs#{'init':'epsg:'+str(epsg)}
        a['area'] = [block[1] for block in outblocks]
        a['perim'] = [block[2] for block in outblocks]
        a['lemgth'] = [block[3] for block in outblocks]
        a['density'] = [1000000/block[1] for block in outblocks]
        b = a.to_crs(epsg=4326)
        b.to_file(folder_name+'blocks'+'latlon'+'.geojson', driver='GeoJSON')
        #b.to_file(folder_name+'blocks'+'latlon'+'.shp')
        filtered_blocks = []
        for block in outblocks:
            if block[3] < 50:
                if 1000 < block[1] < 1000000:
                    filtered_blocks.append(block)
        c = gpd.GeoDataFrame(geometry=[block[0] for block in outblocks])
        c.crs = crs#{'init':'epsg:'+str(epsg)}
        c['area'] = [block[1] for block in outblocks]
        c['perim'] = [block[2] for block in outblocks]
        c['lemgth'] = [block[3] for block in outblocks]
        try:
            blockmedian = statistics.median([block[1] for block in filtered_blocks])
            print('median block density')
            results['blockmedian_density'] = 1000000 / blockmedian
            print(results['blockmedian_density'])
        except:
            results['blockmedian_density'] = 0
            print('BLOCK MEDIAN COULD NOT BE CALCULATED')
        try:
            blockmean = statistics.mean([block[1] for block in filtered_blocks])
            print('mean block density')
            results['blockmean_density'] = 1000000 / blockmean
            print(results['blockmean_density'])
        except:
            results['blockmean_density'] = 0
            print('BLOCK MEAN COULD NOT BE CALCULATED')
        
    ft = datetime.datetime.now()
    print("total", str(ft-dt))
    results['calctime'] = str(ft-dt)
    results['total_pop'] = total_pop
    with open(folder_name+"results.json","w") as output:
        output.write(json.dumps(results))
    for file in ['city.o5m','cityhighways.o5m','citywalk.o5m']:
        if os.path.exists(folder_name+file):
            os.remove(folder_name+file)
    return results

