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

import osmnx as ox
import networkx as nx
import geopandas as gpd
import shapely.geometry
import shapely.ops

import isochrones
import get_service_locations
import car_free_streets
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
    
    patch_length = 5 #kilometers
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
    
    return patches

def weighted_pop_density(array):
    total = 0
    for cell in array.compressed():
        total += cell**2
    return total / numpy.sum(array)

def pedestrians_first(city, 
                      folder_name='', 
                      buffer_dist=100,#m
                      headway_threshold=10,#min
                      to_test = [
                           #'healthcare',
                           #'schools',
                           #'h+s',
                           #'libraries',
                           #'carfree',
                           'blocks',
                           #'density',
                           #'transit',
                           ],
                      distances = { #network buffers, in meters
                            'healthcare': 1000,
                            'schools': 1000,
                            'libraries': 1000,
                            'transit': 500,
                            },
                      overpass = False,
                      patch_length = 5, #km
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
    
    
    boundaries = shapely.geometry.shape(city['geometry'])
    total_pop = city['properties']['P15']
    name = city['properties']['UC_NM_MN']
    hdc = city['properties']['ID_HDC_G0']
    bbox = (city['properties']['BBX_LATMN'],
               city['properties']['BBX_LONMN'],
               city['properties']['BBX_LATMX'],
               city['properties']['BBX_LONMX'],)
    
    
    crs = None 
    
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
    
    patch_times = []
    
    all_coords={}
    
    testing_services = []
    for service in ['healthcare', 'schools', 'libraries']:
        if service in to_test:
            testing_services.append(service)
            #all_coords[service] = get_point_locations(boundaries, queries[service])
    
    if len(testing_services) > 0:
        handler = get_service_locations.ServiceHandler()
        handler.apply_file(str(hdc)+'/city.o5m', locations=True)
        for service in testing_services:
            all_coords[service] = handler.locationlist[service]
            citywide_carfree = handler.carfreelist

    if 'transit' in to_test:
        testing_services.append('transit')
        sources = gtfs_parser.get_feed_infos(gtfs_parser.get_relevant_locs(bbox))
        transit_stop_sets = gtfs_parser.count_all_sources(sources, headwaylim = headway_threshold * 2)
    
    if len(to_test) > 0 and to_test != ["blocks"]:
        for p_idx, patch in enumerate(patches):
            try:
                patch_start = datetime.datetime.now()
                
                unbuffered_patch = patch
                
                max_service_dist_km = max(distances.values())/1000
                
                patch = shapely.geometry.box(
                        patch.bounds[0] - (max_service_dist_km * longitude_factor),
                        patch.bounds[1] - (max_service_dist_km * latitude_factor),
                        patch.bounds[2] + (max_service_dist_km * longitude_factor),
                        patch.bounds[3] + (max_service_dist_km * latitude_factor)
                        )
                        
                
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
                    boundingarg = '-b='
                    boundingarg += str(patch.bounds[0])+','
                    boundingarg += str(patch.bounds[1])+','
                    boundingarg += str(patch.bounds[2])+','
                    boundingarg += str(patch.bounds[3])
                    subprocess.check_call(['osmconvert',
                                           str(hdc)+'/citywalk.o5m',
                                           boundingarg,
                                           #'--complete-ways',
                                           '--drop-broken-refs',
                                           '-o=patch.osm'])
                    G = ox.graph_from_file('patch.osm', 
                                           simplify=False, retain_all=True)
                    os.remove('patch.osm')
                
                G.remove_nodes_from(list(nx.isolates(G)))
                
                simple_G = ox.simplify_graph(G)
                center_nodes = {}
                for service in all_coords.keys():
                    if service in ['healthcare','schools','libraries']:
                        already_covered = None
                        center_nodes[service] = []
                        for coord in all_coords[service]:
                            lat = coord[0]
                            lon = coord[1]
                            if unbuffered_patch.bounds[0] < lon < unbuffered_patch.bounds[2] and unbuffered_patch.bounds[1] < lat < unbuffered_patch.bounds[3]:
                                point = shapely.geometry.Point(lon,lat)
                                if not already_covered:
                                    nearest = ox.get_nearest_node(simple_G, coord)
                                    if not nearest in center_nodes[service]:    
                                        center_nodes[service].append(nearest)
                                    already_covered = point.buffer(50*longitude_factor_m)
                                elif not already_covered.contains(point):
                                    nearest = ox.get_nearest_node(simple_G, coord)
                                    if not nearest in center_nodes[service]:    
                                        center_nodes[service].append(nearest)
                                        already_covered = already_covered.union(point.buffer(50*longitude_factor_m))
                                else:
                                    print ("Already covered")
                
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
                
                
                # Get polygons
                for service in testing_services:
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
        
            except ox.EmptyOverpassResponse:
                print("RECEIVED NO NETWORK FOR PATCH", p_idx)
    
    
                
    epsg = 32600+int(G.graph['crs'].split(' ')[1].split('=')[1]) #This is wild -- 
    #osmnx seems to just give all data in northern-hemisphere format
    #Sorry about the stupid parsing of the projection definition, I'm lazy 
    
    results = {'name':name,'hdc':hdc}
    
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
     
    boundaries_latlon = gpd.GeoDataFrame(geometry=[boundaries])
    boundaries_latlon.crs = {'init':'epsg:4326'}
    boundaries_utm = boundaries_latlon.to_crs(crs)
    
    for service in testing_services:
        if quilt_isochrone_polys[service]:
            service_utm = gpd.GeoDataFrame(geometry = [quilt_isochrone_polys[service]])
            service_utm.crs = {'init':'epsg:'+str(epsg)}
            service_utm.geometry = service_utm.geometry.simplify(15) #maybe this should be after the population calculation
            service_utm = gpd.overlay(service_utm ,boundaries_utm, how='intersection')
            service_latlon = service_utm.to_crs(epsg=4326)
            service_latlon.to_file(folder_name+service+'latlon'+'.geojson', driver='GeoJSON')
            
            stats = rasterstats.zonal_stats(service_latlon, 'pop_dens.tif', stats=['mean'])
            
            total_PNS = 0
            for i in range(0,len(stats)):
                if stats[i]['mean'] and type(stats[i]['mean']) != numpy.ma.core.MaskedConstant:
                    total_PNS += (service_utm.area[i]*stats[i]['mean'] / 62500) #62500 = m2 per pixel
            print("\n")
            print('Total People Near Service for', service, ":", total_PNS)
            print(100*total_PNS/total_pop,"% of",total_pop)
            results[service] = total_PNS / total_pop
        else:
            print ('NO SERVICE FOR', service)
            results[service] = 0
    
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
            
            stats = rasterstats.zonal_stats(carfree_latlon, 'pop_dens.tif', stats=['mean'])
            total_carfree = 0
            for i in range(0,len(stats)):
                if stats[i]['mean'] and type(stats[i]['mean']) != numpy.ma.core.MaskedConstant:
                    total_carfree += (carfree_utm.area[i]*stats[i]['mean'] / 62500) #62500 = m2 per pixel
            print("\n")
            print('Total People Near Service for carfree', ":", total_carfree)
            print(100*total_carfree/total_pop,"% of",total_pop)
            results['carfree'] = total_carfree / total_pop
            
            carfree_latlon = carfree_utm.to_crs('epsg:4326')
            carfree_latlon.to_file(folder_name+'carfreelatlon'+'.geojson', driver='GeoJSON')
            
    
    if 'h+s' in to_test:
        if quilt_isochrone_polys['healthcare'] and quilt_isochrone_polys['schools']:
            service = 'h+s'
            intersect = quilt_isochrone_polys['healthcare'].intersection(quilt_isochrone_polys['schools'])
            if type(intersect) == shapely.geometry.collection.GeometryCollection:
                intersect = [obj for obj in intersect if type(obj) == shapely.geometry.polygon.Polygon]
                intersect = shapely.geometry.MultiPolygon(intersect)
            hs_utm = gpd.GeoDataFrame(geometry = [intersect])
            hs_utm.crs = {'init':'epsg:'+str(epsg)}
            carfree_utm = gpd.overlay(carfree_utm ,boundaries_utm, how='intersection')
            hs_utm.geometry = hs_utm.geometry.simplify(15) #maybe this should be after the population calculation
            hs_latlon = hs_utm.to_crs(epsg=4326)
            hs_latlon.to_file(folder_name+service+'latlon'+'.geojson', driver='GeoJSON')
            
            stats = rasterstats.zonal_stats(hs_latlon, 'pop_dens.tif', stats=['mean'])
            
            total_PNS = 0
            for i in range(0,len(stats)):
                if stats[i]['mean'] and type(stats[i]['mean']) != numpy.ma.core.MaskedConstant:
                    total_PNS += (hs_utm.area[i]*stats[i]['mean'] / 62500) #62500 = m2 per pixel
            print("\n")
            print('Total People Near Service for', service, ":", total_PNS)
            print(100*total_PNS/total_pop,"% of",total_pop)
            results[service] = total_PNS / total_pop
    
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
        
        pdb.set_trace()
        
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
                                       str(hdc)+'/citywalk.o5m',
                                       boundingarg,
                                       #'--complete-ways',
                                       '--drop-broken-refs',
                                       '-o=patch.osm'])
                try:
                    G = ox.graph_from_file('patch.osm', simplify=False, retain_all=True)
                except ox.core.EmptyOverpassResponse:
                    G = False
            
            if G:
                G = ox.project_graph(G, to_crs=crs)
                G = ox.simplify_graph(G)
                
                streets = ox.save_load.graph_to_gdfs(G, nodes = False)
                
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
                                    block = block.simplify(15)
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
            else:
                block_counts.append(0)
        
        #export            
        
        patch_densities = gpd.GeoDataFrame(geometry = list(patches))
        patch_densities['block_count'] = block_counts
        patch_densities.crs = {'init':'epsg:4326'}
        patch_densities_utm = patch_densities.to_crs(epsg=epsg)
        patch_densities_utm['density'] = patch_densities_utm.block_count / (patch_densities_utm.area /1000000)
        patch_densities_latlon = patch_densities_utm.to_crs(epsg=4326)
        patch_densities_latlon.to_file(folder_name+'patch_densities'+'latlon'+'.geojson', driver='GeoJSON')
        
        a = gpd.GeoDataFrame(geometry=[block[0] for block in outblocks])
        a.crs = {'init':'epsg:'+str(epsg)}
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
        c.crs = {'init':'epsg:'+str(epsg)}
        c['area'] = [block[1] for block in outblocks]
        c['perim'] = [block[2] for block in outblocks]
        c['lemgth'] = [block[3] for block in outblocks]
        blockmedian = statistics.median([block[1] for block in filtered_blocks])
        print('median block density')
        results['blockmedian_density'] = 1000000 / blockmedian
        print(results['blockmedian_density'])
        blockmean = statistics.mean([block[1] for block in filtered_blocks])
        print('mean block density')
        results['blockmean_density'] = 1000000 / blockmean
        print(results['blockmean_density'])
        
    ft = datetime.datetime.now()
    print("total", str(ft-dt))
    results['calctime'] = str(ft-dt)
    results['total_pop'] = total_pop
    with open(folder_name+"results.json","w") as output:
        output.write(json.dumps(results))
    for file in ['city.o5m','cityhighways.o5m','citywalk.o5m']:
        if os.path.exists(str(hdc)+'/'+file):
            os.remove(str(hdc)+'/'+file)
    return results

