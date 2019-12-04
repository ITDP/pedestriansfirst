import fiona
import datetime
import os
import numpy
import math
import statistics
import rasterstats
import rasterio.mask
import subprocess

import osmnx as ox
import geopandas as gpd
import networkx as nx
import shapely.geometry
import shapely.ops

import local_isometric
from get_service_locations import queries, get_point_locations
import car_free_streets
import gtfs_parser

import pdb

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
    
    latitude_factor = 0.00898
    longitude_factor = 1/(math.cos(0.0174533 * boundaries.bounds[1]))
    
    height_degrees = abs(boundaries.bounds[3]-boundaries.bounds[1])
    height_km = height_degrees / latitude_factor
    width_degrees = abs(boundaries.bounds[2]-boundaries.bounds[0])
    width_km = width_degrees / longitude_factor
    
    patch_length = 2 #kilometers
    n_vslicers = math.floor(height_km / patch_length)
    n_vslicers = math.floor(width_km / patch_length)
    
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
            all_coords[service] = get_point_locations(boundaries, queries[service])
    
    if 'transit' in to_test:
        testing_services.append('transit')
        sources = gtfs_parser.get_feed_infos(gtfs_parser.get_relevant_locs(bbox))
        transit_stop_sets = gtfs_parser.count_all_sources(sources, headwaylim = headway_threshold * 2)
        
    if 'blocks' in to_test:
        outblocks = [] 
        
    for p_idx, patch in enumerate(patches):
        try:
            patch_start = datetime.datetime.now()
            
            if abs(patch.bounds[1]) < 23:
                patch_buffer = .0055
            elif abs(patch.bounds[1]) < 45:
                patch_buffer = .0067
            elif abs(patch.bounds[1]) < 67:
                patch_buffer = .0114
            else:
                patch_buffer = .02 #in decimal degrees
            
            unbuffered_patch = patch
            patch = patch.buffer(patch_buffer)
            
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
            for service in all_coords.keys():
                if service in ['healthcare','schools','libraries']:
                    center_nodes[service] = []
                    for coord in all_coords[service]:
                        lat = coord[0]
                        lon = coord[1]
                        if patch.bounds[0] < lon < patch.bounds[2] and patch.bounds[1] < lat < patch.bounds[3]:
                            nearest = ox.get_nearest_node(simple_G, coord)
                            if not nearest in center_nodes[service]:
                                center_nodes[service].append(nearest)
            
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
            
            #we have a crs now, so:
            unbuffered_patch = gpd.GeoSeries(unbuffered_patch, crs={'init':'epsg:4326'}).to_crs(crs)[0]
            
            isochrone_polys = {}
            failures = {}
            for service in to_test:
                failures[service] = 0
        
            if 'carfree' in to_test:
                polygons = []
                ped_G, park_G = car_free_streets.get_carfree_graph(patch)
                if ped_G and park_G:
                    carfree_G = ox.project_graph(nx.compose(ped_G, park_G), to_crs = crs)
                if ped_G:
                    carfree_G = ox.project_graph(ped_G, to_crs = crs)
                if park_G:
                    carfree_G = ox.project_graph(park_G, to_crs = crs)
                else:
                    carfree_G = None
                if carfree_G:
                    try:
                        isochrone_polys['carfree'] = local_isometric.network_buffer(carfree_G, distance=buffer_dist)
                    except AttributeError:
                        isochrone_polys['carfree'] = None
            
            if 'blocks' in to_test:
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
                                    filtered_blocks.append(block)
                        outblocks += filtered_blocks
                    
                
                
            
            # Get polygons
            for service in testing_services:
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
    
    #Export    
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
    
    #pdb.set_trace()
    
    if 'blocks' in to_test:
        a = gpd.GeoDataFrame(geometry=outblocks)
        a.crs = {'init':'epsg:'+str(epsg)}
        area = a.geometry.area
        perim = a.geometry.length
        a['area'] = area
        a['perim'] = perim
        a['lemgth'] = (area*area)/perim
        a.geometry = a.geometry.simplify(15)
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
