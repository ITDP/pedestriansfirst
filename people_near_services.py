import fiona
import datetime
import os
import numpy
import math
import statistics
import rasterstats
import rasterio.mask
import subprocess
import osmium

import osmnx as ox
import geopandas as gpd
import networkx as nx
import shapely.geometry
import shapely.ops
import pyproj

import local_isometric
import get_service_locations
import car_free_streets
import gtfs_parser

import pdb
import logging


ox.utils.config(log_console = False)

#added a test comment and changed it

def weighted_pop_density(array):
    total = 0
    for cell in array.compressed():
        total += cell**2
    return total / numpy.sum(array)

def pnservices(city, folder_name='', buffer_dist=100, headway_threshold=10,
               to_test = [
                       'healthcare',
                       'schools',
                       'h+s',
                       'libraries',
                       'carfree',
                       'blocks',
                       'density',
                       'transit',
                       ],
                distances = {
                        'healthcare': 1000,
                        'schools': 1000,
                        'libraries': 1000,
                        'transit': 500,
                        },
                overpass = False
                ):    
    dt = datetime.datetime.now()

#shpfile_location = 'rva.shp'
#folder_name = 'rva\\'
#network_dist = 500
#buffer_dist = 100
#headway_threshold = 20
#to_test = [#'healthcare',
#           #'schools',
#           #'libraries',
#           #'carfree',
#           #'parks'
#           'blocks',
#           #'density',
#           #'transit',
#           #'transit2',
#           ]
    dt = datetime.datetime.now()
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
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
    
    SAVED_TSTOPS = []
    
    crs = None    
    
    print('Evaluating Pedestrians First indicators in',name)
    print('Measuring',str(to_test))
    
    
    longitude_factor = 0.00898 # degrees per km
    longitude_factor_m = 0.00898 / 1000 # degrees per m
    latitude_factor = 1/(math.cos(boundaries.bounds[1]/0.0174533)) 
    
    height_degrees = abs(boundaries.bounds[3]-boundaries.bounds[1])
    height_km = height_degrees / latitude_factor
    width_degrees = abs(boundaries.bounds[2]-boundaries.bounds[0])
    width_km = width_degrees / longitude_factor
    
    patch_length = 2 #kilometers
    n_vslicers = math.floor(height_km / patch_length)
    n_hslicers = math.floor(width_km / patch_length)
    
    hslicers = []
    vslicers = []
    
    for i in range(1,n_hslicers+1):
        increment = (bbox[2]-bbox[0])/(n_hslicers+1)
        lat = bbox[0]+(i*increment)
        slicer = shapely.geometry.LineString([(bbox[1],lat),(bbox[3],lat)])
        hslicers.append(slicer)
        
    for i in range(1,n_vslicers+1):
        increment = (bbox[3]-bbox[1])/(n_vslicers+1)
        lon = bbox[1]+(i*increment)
        slicer = shapely.geometry.LineString([(lon,bbox[0]),(lon, bbox[2])])
        hslicers.append(slicer)
        
    patches = shapely.geometry.MultiPolygon(polygons=[boundaries])
    for slicer in hslicers+vslicers:
        patches = shapely.geometry.MultiPolygon(polygons = shapely.ops.split(patches, slicer))
    
    
    
    
    for i in range(1,n_hslicers+1):
        h_increment = (bbox[2]-bbox[0])/(n_hslicers+1)
        lat = bbox[0]+(i*h_increment)
        slicer = shapely.geometry.LineString([(bbox[1],lat),(bbox[3],lat)])
        hslicers.append(slicer)
        
    for i in range(1,n_vslicers+1):
        v_increment = (bbox[3]-bbox[1])/(n_vslicers+1)
        lon = bbox[1]+(i*v_increment)
        slicer = shapely.geometry.LineString([(lon,bbox[0]),(lon, bbox[2])])
        hslicers.append(slicer)
        
    patches = shapely.geometry.MultiPolygon(polygons=[boundaries])
    for slicer in hslicers+vslicers:
        patches = shapely.geometry.MultiPolygon(polygons = shapely.ops.split(patches, slicer))
    
    print("cut"+str(len(list(patches)))+"patches")
    
    quilt_ipolys = {}
    quilt_cnodes = {}
    for service in to_test:
        quilt_ipolys[service] = False
        quilt_cnodes[service] = []
    
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
                    
            
            walk_filter = ('["area"!~"yes"]["highway"!~"link|motor|proposed|construction|abandoned|platform|raceway"]'
                           '["service"!~"parking_aisle|driveway"]'
                                   '["foot"!~"no"]["service"!~"private"]{}').format(ox.settings.default_access)
            
            if overpass:
                G = ox.graph_from_polygon(patch, custom_filter=walk_filter, simplify=False)
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
                G = ox.graph_from_file('patch.osm', simplify=False)
                #os.remove('patch.osm')
            
            simple_G = ox.simplify_graph(G)
            center_nodes = {}
            print("getting service center_node")
            for service in all_coords.keys():
                if service in ['healthcare','schools','libraries']:
                    already_covered = None
                    center_nodes[service] = []
                    for coord in all_coords[service]:
                        lat = coord[0]
                        lon = coord[1]
                        if patch.bounds[0] < lon < patch.bounds[2] and patch.bounds[1] < lat < patch.bounds[3]:
                            point = shapely.geometry.Point(lon,lat)
                            if not already_covered:
                                already_covered = point.buffer(50*longitude_factor_m)
                            elif not already_covered.contains(point):
                                nearest = ox.get_nearest_node(simple_G, coord)
                                if not nearest in center_nodes[service]:    
                                    center_nodes[service].append(nearest)
                                    already_covered = already_covered.union(point.buffer(50*longitude_factor_m))
            
            print('getting transit center_nodes')
            if 'transit' in to_test:
                center_nodes['transit'] = []
                transit_centers = {}
                for service_idx, service in enumerate(transit_stop_sets):
                    for stop_id in service.keys():
                        lat = float(service[stop_id][1])
                        lon = float(service[stop_id][2])
                        headway = float(service[stop_id][0])
                        if patch.bounds[0] < lon < patch.bounds[2] and patch.bounds[1] < lat < patch.bounds[3]:
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
                            if min_hw > headway:
                                SAVED_TSTOPS.append((min_hw, headway))
                            
            
            
            
            
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
            
            print("getting carfree")
            if 'carfree' in to_test:
                carfree = []
                for poly in citywide_carfree:
                    rep = poly.representative_point() 
                    if patch.bounds[0] < rep.x < patch.bounds[2] and patch.bounds[1] < rep.y < patch.bounds[3]:
                        carfree.append(poly)
                carfree = shapely.ops.cascaded_union(carfree)
                if carfree:
                    pdb.set_trace()
                    print(crs)
                    projection = pyproj.Transformer.from_crs(4326, crs)
                    carfree = shapely.ops.transform(projection.transform, carfree)
                    isochrone_polys['carfree'] = carfree.buffer(buffer_dist)
                
            
            
                    
            
            
            # Get polygons
            for service in testing_services:
                print(service)
                isochrone_polys[service], fails = local_isometric.make_iso_polys(G, center_nodes[service], distance=distances[service], edge_buff=buffer_dist)
                failures[service] += fails
                
            for service in isochrone_polys.keys():
                if service not in quilt_ipolys.keys() or not quilt_ipolys[service]:
                    quilt_ipolys[service] = isochrone_polys[service]
                elif isochrone_polys[service]:
                    quilt_ipolys[service] = shapely.ops.cascaded_union([quilt_ipolys[service],isochrone_polys[service]])
            for service in center_nodes.keys():
                quilt_cnodes[service] = quilt_cnodes[service] + center_nodes[service]
            
            patch_time = datetime.datetime.now() - patch_start
                
            print("finished patch #",p_idx,'out of', len(patches),"in",str(patch_time))
            patch_times.append(patch_time)
    
        except ox.EmptyOverpassResponse:
            print("RECEIVED NO NETWORK FOR PATCH", p_idx)
    
    
                
    epsg = 32600+int(G.graph['crs'].split(' ')[1].split('=')[1]) #This is wild -- 
    #osmnx seems to just give all data in northern-hemisphere format
    #Sorry about the stupid parsing of the projection definition, I'm lazy 
    
    results = {}
    
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
       
    if 'carfree' in to_test:
        testing_services.append('carfree')
    
    for service in testing_services:
        if quilt_ipolys[service]:
            a = gpd.GeoDataFrame(geometry = [quilt_ipolys[service]])
            a.crs = {'init':'epsg:'+str(epsg)}
            a.geometry = a.geometry.simplify(15) #maybe this should be after the population calculation
            b = a.to_crs(epsg=4326)
            b.to_file(folder_name+service+'latlon'+'.geojson', driver='GeoJSON')
            b.to_file(folder_name+service+'latlon'+'.shp') #unnecessary later
            
            #a, b = local_isometric.export(quilt_ipolys[service], epsg, service=service, folder=folder_name)
            
            stats = rasterstats.zonal_stats(b, 'pop_dens.tif', stats=['mean'])
            
            total_PNS = 0
            for i in range(0,len(stats)):
                if stats[i]['mean'] and type(stats[i]['mean']) != numpy.ma.core.MaskedConstant:
                    total_PNS += (a.area[i]*stats[i]['mean'] / 62500) #62500 = m2 per pixel
            print("\n")
            print('Total People Near Service for', service, ":", total_PNS)
            print(100*total_PNS/total_pop,"% of",total_pop)
            results[service] = total_PNS / total_pop
        else:
            print ('NO SERVICE FOR', service)
            results[service] = 0
    
    if 'h+s' in to_test:
        if quilt_ipolys['healthcare'] and quilt_ipolys['schools']:
            service = 'h+s'
            intersect = quilt_ipolys['healthcare'].intersection(quilt_ipolys['schools'])
            if type(intersect) == shapely.geometry.collection.GeometryCollection:
                intersect = [obj for obj in intersect if type(obj) == shapely.geometry.polygon.Polygon]
                intersect = shapely.geometry.MultiPolygon(intersect)
            a = gpd.GeoDataFrame(geometry = [intersect])
            a.crs = {'init':'epsg:'+str(epsg)}
            a.geometry = a.geometry.simplify(15) #maybe this should be after the population calculation
            b = a.to_crs(epsg=4326)
            b.to_file(folder_name+service+'latlon'+'.geojson', driver='GeoJSON')
            try:
                b.to_file(folder_name+service+'latlon'+'.shp') #unnecessary later
            except:
                pdb.set_trace()
                
            #a, b = local_isometric.export(quilt_ipolys[service], epsg, service=service, folder=folder_name)
            
            stats = rasterstats.zonal_stats(b, 'pop_dens.tif', stats=['mean'])
            
            total_PNS = 0
            for i in range(0,len(stats)):
                if stats[i]['mean'] and type(stats[i]['mean']) != numpy.ma.core.MaskedConstant:
                    total_PNS += (a.area[i]*stats[i]['mean'] / 62500) #62500 = m2 per pixel
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
    quilt_ipolys = None
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
        #most of our patch-cutting variables are still around
        patch_length = .5 #kilometers
        n_vslicers = math.floor(height_km / patch_length)
        n_hslicers = math.floor(width_km / patch_length)
        
        hslicers = []
        vslicers = []
        
        for i in range(1,n_hslicers+1):
            increment = (bbox[2]-bbox[0])/(n_hslicers+1)
            lat = bbox[0]+(i*increment)
            slicer = shapely.geometry.LineString([(bbox[1],lat),(bbox[3],lat)])
            hslicers.append(slicer)
            
        for i in range(1,n_vslicers+1):
            increment = (bbox[3]-bbox[1])/(n_vslicers+1)
            lon = bbox[1]+(i*increment)
            slicer = shapely.geometry.LineString([(lon,bbox[0]),(lon, bbox[2])])
            hslicers.append(slicer)
            
        patches = shapely.geometry.MultiPolygon(polygons=[boundaries])
        for slicer in hslicers+vslicers:
            patches = shapely.geometry.MultiPolygon(polygons = shapely.ops.split(patches, slicer))
        
        print ("cut", len(list(patches)),"patches")
        n=0
        outblocks = []
        for patch in patches:
            print("patch"+str(n)+" of "+str(len(patches)) )
            n+=1
            unbuffered_patch = gpd.GeoSeries(patch, crs={'init':'epsg:4326'})
            unbuffered_patch = unbuffered_patch.to_crs(crs)[0]
            
            
            buffer_dist = .2
            
            patch = shapely.geometry.box(
                    patch.bounds[0] - (buffer_dist * longitude_factor),
                    patch.bounds[1] - (buffer_dist * latitude_factor),
                    patch.bounds[2] + (buffer_dist * longitude_factor),
                    patch.bounds[3] + (buffer_dist * latitude_factor)
                    )
            
            if overpass:
                G = ox.graph_from_polygon(patch, custom_filter=walk_filter, simplify=False)
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
                G = ox.graph_from_file('patch.osm', simplify=False)
                    
            G = ox.project_graph(G, to_crs=crs)
            G = ox.simplify_graph(G)
            
            streets = ox.save_load.graph_to_gdfs(G, nodes = False)
            
            if not streets.empty:
                streets = shapely.geometry.MultiLineString(list(streets.geometry))
                merged = shapely.ops.linemerge(streets)
                if merged:
                    borders = shapely.ops.unary_union(merged)
                    blocks = list(shapely.ops.polygonize(borders))
                    filtered_blocks = []
                    for block in blocks:
                        if 1000 < block.area < 1000000:
                            if block.interiors:
                                block = shapely.geometry.Polygon(block.exterior)
                            #if block.length / block.area < 0.15:
                            if block.centroid.within(unbuffered_patch):
                                area = round(block.area, 3)
                                perim = round(block.length, 3)
                                lemgth = round((perim * perim) / area, 3)
                                if lemgth < 50:
                                    block = block.simplify(15)
                                    filtered_blocks.append((block, area, perim, lemgth))
                    outblocks += filtered_blocks  
        
        #export            
        
        a = gpd.GeoDataFrame(geometry=[block[0] for block in outblocks])
        a.crs = {'init':'epsg:'+str(epsg)}
        a['area'] = [block[1] for block in outblocks]
        a['perim'] = [block[2] for block in outblocks]
        a['lemgth'] = [block[3] for block in outblocks]
        b = a.to_crs(epsg=4326)
        b.to_file(folder_name+'blocks'+'latlon'+'.geojson', driver='GeoJSON')
        b.to_file(folder_name+'blocks'+'latlon'+'.shp')
        blockmedian = statistics.median(a.area)
        print('median block size')
        print(blockmedian)
        results['blockmedian'] = blockmedian
        blockmean = statistics.mean(a.area)
        print('mean block size')
        print(blockmean)
        results['blockmean'] = blockmean
        
    
    
    ft = datetime.datetime.now()
    print("total", str(ft-dt))
    with open(folder_name+"results.txt","w") as output:
        for item in results.keys():    
            output.write(item+": "+str(results[item])+"\n")
    return results
#return failures, SAVED_TSTOPS

#failures, saved_tstops = pnservices('dc.shp',folder_name='dc\\')
#CITIES = ["atlanta","nyc","houston","mexico","indianapolis","montreal","guadalajara","leon","monterrey","washington","toronto","ottowa","atlanta"]
#x = pnservices('arax.shp',folder_name='arax/')
