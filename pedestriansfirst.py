import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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
import shutil
#import utm_zone

import numpy as np
import osmnx as ox
import networkx as nx
import pandas as pd
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
import geopandas as gpd
import shapely.geometry
from shapely.geometry import LineString, Point
from shapely.ops import unary_union
import shapely.ops

import isochrones
import get_service_locations
import gtfs_parser

import pdb


ox.utils.config(log_console = False)

def make_patches(bound_latlon, crs_utm, patch_length = 10000, buffer = 500): #patch_length and buffer in m
    ''' 
    'tile' the boundaries of a city into patches, like a patchwork quilt, including a buffer
    '''
    bounds_utm_poly = bound_latlon.to_crs(crs_utm).geometry.unary_union
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
    
    print(f"cut {len(patches)} patches")
    
    return patches_latlon

def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
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

def add_intermediate_nodes(G, maxdist):
    start_time = datetime.datetime.now()
    print(f'adding nodes to a max edge length of {maxdist}')
    nodes, edges = ox.graph_to_gdfs(G)
    newnode_baseid = nodes.index.max()
    nodes_added = 0
    # tags_to_copy = list(edges.columns)
    # tags_to_copy.remove('length')
    # tags_to_copy.remove('geometry')
    while len(edges[edges.geometry.length>maxdist]) > 0:
        to_shorten = edges[edges.geometry.length>maxdist].iloc[0]
        
        old_leg = to_shorten.geometry
        
        new_legs, new_node = cut(old_leg, old_leg.length/2)
        
        new_node_id = newnode_baseid + nodes_added
        nodes_added += 1
        nodes.loc[new_node_id, 'y'] = new_node.y
        nodes.loc[new_node_id, 'x'] = new_node.x
        nodes.loc[new_node_id, 'geometry'] = new_node
        
        first_leg_id = (to_shorten.name[0], new_node_id, to_shorten.name[2])
        edges.loc[first_leg_id,'geometry'] = new_legs[0]
        edges.loc[first_leg_id,'length'] = new_legs[0].length
        second_leg_id = (new_node_id, to_shorten.name[1], to_shorten.name[2])
        edges.loc[second_leg_id,'geometry'] = new_legs[1]
        edges.loc[second_leg_id,'length'] = new_legs[1].length
        
        # for tag in tags_to_copy:
        #     edges.loc[first_leg_id, tag] = to_shorten.loc[tag]
        #     edges.loc[second_leg_id, tag] = to_shorten.loc[tag]
        
        edges = edges.drop(to_shorten.name)
    end_time = datetime.datetime.now()
    print(f"added {nodes_added} nodes in {end_time - start_time}")
    return ox.graph_from_gdfs(nodes, edges)
        

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
                           #'bikeshare',
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
                      patch_length = 10000, #m
                      boundary_buffer = 500, #m
                      blocks_simplification = 15, #m
                      gtfs_files = [],
                      debug = False,
                      ):   
    
    dt = datetime.datetime.now()
    
    
    useful_tags = ox.settings.useful_tags_way + ['cycleway', 'cycleway:left', 'cycleway:right', 'cycleway:both', 'bicycle']
    ox.config(use_cache=True, log_console=True, useful_tags_way=useful_tags)
    
    if folder_name != '' and not folder_name[-1:] == '/':
        folder_name += '/'
    
    bound_latlon = gpd.GeoDataFrame(geometry = [boundaries])
    bound_latlon.crs = {'init':'epsg:4326'}
    longitude = round(numpy.mean(bound_latlon.geometry.centroid.x),10)
    utm_zone = int(math.floor((longitude + 180) / 6) + 1)
    utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)
    
    
    if boundary_buffer > 0:
        bound_utm = bound_latlon.to_crs(utm_crs)
        bound_utm.geometry = bound_utm.geometry.buffer(boundary_buffer)
        bound_latlon = bound_utm.to_crs(epsg=4326)
        boundaries = bound_latlon.geometry.unary_union
    print(utm_crs)
    bbox = boundaries.bounds
    
    
    patches = make_patches(bound_latlon, utm_crs, patch_length=patch_length)
    
    print('Evaluating Pedestrians First indicators in',name)
    print('Measuring',str(to_test))
    
    quilt_isochrone_polys = {}
    for service in to_test:
        quilt_isochrone_polys[service] = False
        if 'pnab' in to_test or 'pnpb' in to_test:
            quilt_isochrone_polys['pnab'] = False
            quilt_isochrone_polys['pnpb'] = False
            
    
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
        
    if len(to_test) > 0 and to_test != ["blocks"]:
        for p_idx, patch in enumerate(patches.geometry):
            try:
                patch_start = datetime.datetime.now()
                
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
                        G = ox.graph_from_xml('patch.osm', simplify=False, retain_all=False)
                        os.remove('patch.osm')
                    except TypeError: #something to do with clipping, seems to happen once in a while
                        #pdb.set_trace()
                        #this is a very stupid band-aid, but it works for now
                        print ('TYPEERROR FROM CLIPPING PATCH', p_idx)
                        with open(str(folder_name)+"patcherrorlog.txt", "a") as patcherrorlog:
                            patcherrorlog.write('TYPEERROR FROM CLIPPING PATCH '+str(p_idx))
                        G = ox.graph_from_polygon(patch, 
                                              custom_filter=walk_filter, 
                                              simplify=False, 
                                              retain_all=True)
                os.remove('patchbounds.poly')
                        
                
                G.remove_nodes_from(list(nx.isolates(G)))
                
                simple_G = ox.simplify_graph(G)
                
                center_nodes = {}
                for service in all_coords.keys():
                    if service in ['healthcare','schools','libraries','bikeshare','special']:
                        center_nodes[service] = []
                        for coord in all_coords[service]:
                            lat = coord[0]
                            lon = coord[1]
                            nearest = ox.get_nearest_node(simple_G, coord)
                            if not nearest in center_nodes[service]:    
                                center_nodes[service].append(nearest)
                
                if 'pnab' in to_test or 'pnpb' in to_test:
                    ways_gdf = ox.graph_to_gdfs(G, nodes=False)
                    
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
                    
                    center_nodes['pnpb'] = set()
                    center_nodes['pnab'] = set()
                    for edge in total_protectedbike.index:
                        if edge[0] in simple_G.nodes:
                            center_nodes['pnpb'].add(edge[0])
                        if edge[1] in simple_G.nodes:
                            center_nodes['pnpb'].add(edge[1])
                    for edge in total_allbike.index:
                        if edge[0] in simple_G.nodes:
                            center_nodes['pnab'].add(edge[0])
                        if edge[1] in simple_G.nodes:
                            center_nodes['pnab'].add(edge[1])
                           
                if 'transit' in to_test:
                    center_nodes['transit'] = []
                    transit_centers = {}
                    for service_idx, service in enumerate(transit_stop_sets):
                        for stop_id in service.keys():
                            lat = float(service[stop_id][1])
                            lon = float(service[stop_id][2])
                            headway = float(service[stop_id][0])
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
                
                
                G = ox.project_graph(simple_G, to_crs=utm_crs)
                G = add_intermediate_nodes(G, 100)
                
                isochrone_polys = {}
                failures = {}
                for service in to_test:
                    failures[service] = 0
                    if 'pnab' in to_test or 'pnpb' in to_test:
                        failures['pnab']=0
                        failures['pnpb']=0
                
                
                
                # Get polygons
                for service in testing_services:
                    print(f'getting polygons for {service}, {len(center_nodes[service])} center_nodes')
                    isochrone_polys[service], fails = isochrones.make_iso_polys(G,
                                                                                center_nodes[service],
                                                                                distance=distances[service],
                                                                                edge_buff=buffer_dist)
                    failures[service] += fails
                    
                for service in isochrone_polys.keys():
                    if service not in quilt_isochrone_polys.keys() or not quilt_isochrone_polys[service]:
                        quilt_isochrone_polys[service] = isochrone_polys[service]
                    elif isochrone_polys[service]:
                        quilt_isochrone_polys[service] = shapely.ops.unary_union([quilt_isochrone_polys[service],isochrone_polys[service]])
                
                patch_time = datetime.datetime.now() - patch_start
                    
                print("finished patch #",p_idx,'out of', len(patches),"in",str(patch_time))
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
    
    results = {'name':name,'id_code':id_code}
    
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
     
    boundaries_latlon = gpd.GeoDataFrame(geometry=[boundaries])
    boundaries_latlon.crs = {'init':'epsg:4326'}
    boundaries_utm = boundaries_latlon.to_crs(utm_crs)
        
    stats = rasterstats.zonal_stats(boundaries_latlon, 'input_data/pop_dens.tif', stats=['sum'])           
    total_pop = stats[0]['sum'] 
    
    for service in testing_services:
        if quilt_isochrone_polys[service]:
            service_utm = gpd.GeoDataFrame(geometry = [quilt_isochrone_polys[service]])
            service_utm.crs = utm_crs#{'init':'epsg:'+str(epsg)}
            service_utm.geometry = service_utm.geometry.simplify(15) #maybe this should be after the population calculation
            service_utm = gpd.overlay(service_utm ,boundaries_utm, how='intersection')
            service_latlon = service_utm.to_crs(epsg=4326)
            service_latlon.to_file(folder_name+service+'latlon'+'.geojson', driver='GeoJSON')
            
            stats = rasterstats.zonal_stats(service_latlon, 'input_data/pop_dens.tif', stats=['sum'])           
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
        if not total_protectedbike.empty:
            merged_protectedbike = gpd.GeoDataFrame(geometry = [total_protectedbike.geometry.unary_union], crs=4326)
            merged_protectedbike = gpd.clip(merged_protectedbike, boundaries)
            merged_protectedbike.crs = 4326
            merged_protectedbike.to_file(folder_name+'bikeways/protectedbike.geojson',driver='GeoJSON')
            protected_km += sum(merged_protectedbike.to_crs(utm_crs).geometry.length)
        if not total_allbike.empty:
            merged_allbike = gpd.GeoDataFrame(geometry = [total_allbike.geometry.unary_union], crs=4326)
            merged_allbike = gpd.clip(merged_allbike, boundaries)
            merged_allbike.crs = 4326
            merged_allbike.to_file(folder_name+'bikeways/allbike.geojson',driver='GeoJSON')
            unprotected_km += sum(merged_allbike.to_crs(utm_crs).geometry.length)
            
        results['protected_bikeways_km'] = protected_km
        results['all_bikeways_km'] = protected_km + unprotected_km
    
    if 'carfree' in to_test:
        print("getting carfree")
        if citywide_carfree:
            carfree_latlon = gpd.GeoDataFrame(geometry = citywide_carfree)
            #just a latlon list of points
            carfree_latlon.crs = {'init':'epsg:4326'}
            carfree_utm = carfree_latlon.to_crs(utm_crs)
            carfree_utm.geometry = carfree_utm.geometry.buffer(100)
            #this is the analysis, the 100m buffer
            carfree_utm = gpd.GeoDataFrame(geometry = [shapely.ops.unary_union(carfree_utm.geometry)])
            carfree_utm.geometry = carfree_utm.geometry.simplify(10)
            carfree_utm = gpd.overlay(carfree_utm ,boundaries_utm, how='intersection')
            carfree_utm.crs = utm_crs
            carfree_latlon = carfree_utm.to_crs('epsg:4326')
            
            stats = rasterstats.zonal_stats(carfree_latlon, 'input_data/pop_dens.tif', stats=['sum'])
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
            hs_utm = gpd.GeoDataFrame(geometry = [intersect], crs=utm_crs)
            if hs_utm.geometry.area.sum() != 0:
                hs_utm = gpd.overlay(hs_utm ,boundaries_utm, how='intersection')
                hs_utm.geometry = hs_utm.geometry.simplify(15)
                hs_latlon = hs_utm.to_crs(epsg=4326)
                hs_latlon.to_file(folder_name+service+'latlon'+'.geojson', driver='GeoJSON')
                stats = rasterstats.zonal_stats(hs_latlon, 'input_data/pop_dens.tif', stats=['sum'])
                
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
                                'input_data/pop_dens.tif', 
                                stats = [],
                                add_stats={'weighted': weighted_pop_density}
                                )[0]['weighted']
        results['density'] = density / 0.0625 #km^2 / pixel
        print('weighted pop density', results['density'])
        with rasterio.open('input_data/pop_dens.tif') as dataset:
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
        for n, patch in enumerate(patches.geometry):
            print("patch"+str(n)+" of "+str(len(patches)) )
            unbuffered_patch = gpd.GeoSeries(patch, crs={'init':'epsg:4326'})
            unbuffered_patch = unbuffered_patch.to_crs(utm_crs)[0]
            
            
           
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
                    G = ox.project_graph(G, to_crs=utm_crs)
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
        patch_densities_utm = patch_densities.to_crs(utm_crs)
        patch_densities_utm['density'] = patch_densities_utm.block_count / (patch_densities_utm.area /1000000)
        patch_densities_latlon = patch_densities_utm.to_crs(epsg=4326)
        patch_densities_latlon.to_file(folder_name+'patch_densities'+'latlon'+'.geojson', driver='GeoJSON')
        
        a = gpd.GeoDataFrame(geometry=[block[0] for block in outblocks])
        a.crs = utm_crs
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
        c.crs = utm_crs
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

