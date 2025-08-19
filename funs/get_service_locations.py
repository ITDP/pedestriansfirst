import requests 
import json
import fiona
import osmnx as ox
import numpy as np
import geopandas as gpd
import shapely
from shapely.geometry import Point

import osmium
import shapely.wkb
import timeit
import gc

# simple_projected_G = G_allroads
# simple_projected_G = ox.load_graphml('teste_kaue/sp_graph_test.graphml')
# min_length = 1000
# min_lanes = 1

# get_highways(simple_projected_G = ox.load_graphml('teste_kaue/sp_graph_test.graphml'))

def get_highways(nodes, edges,
                 min_length = 1000, #meters
                 min_lanes = 1):
    #get the car-only network and the links of major (divided) roads
    # timewise it's okay here
    start = timeit.default_timer()
    
    end = timeit.default_timer()
    end - start
    
    car_tags = ['motorway','trunk','primary','secondary','tertiary','unclassified',
                'residential','living_street','service','road']
    for tag in car_tags.copy():
        car_tags.append(tag+'_link')
    major_tags = ['motorway','trunk','primary','secondary',
                  'motorway_link','trunk_link','primary_link','secondary_link']
    car_roads = edges[edges.highway.isin(car_tags)]
    multi_car_G = ox.graph_from_gdfs(nodes, car_roads)
    major_roads_utm = edges[(edges.highway.isin(major_tags)) & (edges.oneway == True)].copy()
    #only include major roads with at least 2 lanes per direction
    if 'lanes' in major_roads_utm.columns:
        for idx in major_roads_utm.index:
            lanes = major_roads_utm.loc[idx, 'lanes']
            if type(lanes) == type('string'):
                lanes = lanes.split(';')
            if type(lanes) == type([]):
                all_lanes = []
                for lane_ct in lanes:
                    if type(lane_ct) == type('string'):
                        for lane_ct_2 in lane_ct.split(';'):
                            try:
                                all_lanes.append(float(lane_ct_2))
                            except:
                                pass
                    else:
                        all_lanes.append(lane_ct)
                if len(all_lanes) > 0:
                    lanes = min(all_lanes)
                else:
                    lanes = 3
            lanes = float(lanes)
            if np.isnan(lanes):
                lanes = 3 #if the number of lanes isn't given, we assume it's more than 2 per direction
            if lanes < min_lanes: 
                major_roads_utm.drop(idx, inplace=True)
    if len(major_roads_utm) == 0:
        print('len(major_roads_utm) == 0:')
        return None
            
    #exclude long tunnels
    if 'tunnel' in major_roads_utm.columns:
        major_roads_utm.drop(major_roads_utm[(major_roads_utm.tunnel=='yes')&(major_roads_utm.length>300)].index, inplace=True)
    
    
    start = timeit.default_timer()
    # Identify all the nodes with no more than three neighbors 
    # ie, exclude four-way intersections
    major_nodes = set()
    for idx in major_roads_utm.index:
        major_nodes.add(idx[0])
        major_nodes.add(idx[1])
    grade_separated_nodes = []
    at_grade_nodes = []
    for center_node in major_nodes:
        neighbors = set()
        for edge in list(multi_car_G.out_edges(center_node)):
            for neighbor in edge:
                neighbors.add(neighbor)
        for edge in list(multi_car_G.in_edges(center_node)):
            for neighbor in edge:
                neighbors.add(neighbor)
        if center_node in neighbors:
            neighbors.remove(center_node)
        if len(neighbors) > 3:
            at_grade_nodes.append(center_node)
        else:
            grade_separated_nodes.append(center_node)
    end = timeit.default_timer()
    print('Time prep final1: ', end - start)
            
    # split up the highway lines at 4-way intersections
    
    start = timeit.default_timer()
    separation_breakers = nodes.loc[at_grade_nodes]
    separation_break_poly = separation_breakers.buffer(2).unary_union
    separation_break_gdf_utm = gpd.GeoDataFrame(geometry = [separation_break_poly], crs = separation_breakers.crs)
    end = timeit.default_timer()
    print('Time prep final2: ', end - start)
    
    
    if major_roads_utm.empty:
        return major_roads_utm.to_crs(4326)
    
    
    start = timeit.default_timer()    
    # Split major roads at intersections??
    # VERY SLOW HERE!!!!!!!!
    separated_major_roads = major_roads_utm.overlay(separation_break_gdf_utm, how='difference').geometry
    end = timeit.default_timer()
    print('Time prep final3: ', end - start)
    
    #separated_major_roads = separated_major_roads.explode() #boom
    
    #merge lines based on shared nodes
    # all_lines_and_multilines = list(separated_major_roads.geometry)
    # lines_only = []
    # for item in all_lines_and_multilines:
    #     if type(item).__name__ == "LineString":
    #         lines_only.append(item)
    #     elif type(item).__name__ == "MultiLineString":
    #         for subitem in item.explode():
    #             lines_only.append(subitem)
    
    try:
        merged_lines = shapely.ops.linemerge(list(separated_major_roads.geometry))
    except: #if one if the rows has a multilinestring
        merged_lines = shapely.ops.linemerge(list(separated_major_roads.explode().geometry))
        
    try:
        merged_lines_gdf = gpd.GeoDataFrame(geometry=list(merged_lines.geoms), crs=major_roads_utm.crs)
    except AttributeError: #it's a single LineString
        merged_lines_gdf = gpd.GeoDataFrame(geometry=[merged_lines], crs=major_roads_utm.crs)
        
    
    #see which line segments are part of a touching network with a total length greater than the minimum
    start = timeit.default_timer()    
    real_highway_lineids = set()
    for idx in merged_lines_gdf.index:
        #is it already accounted for?
        if idx in real_highway_lineids:
            continue
        #is it, alone, long enough?
        this_segment_termini = merged_lines_gdf.loc[idx,'geometry'].boundary.geoms
        if len(this_segment_termini) < 2: #it's a circle
            if merged_lines_gdf.loc[idx,'geometry'].length > min_length:
                real_highway_lineids.update([idx]) 
            continue
        if this_segment_termini[0].distance(this_segment_termini[1]) > min_length:
            real_highway_lineids.update([idx]) 
            continue
        #or does it touch other segments that are cumulatively long enough?
        contiguous_ids = [idx]
        contiguous_termini = []
        touches_a_line_thats_already_included = False
        for connecting_idx in contiguous_ids:
            immediately_touching = list(merged_lines_gdf[merged_lines_gdf.geometry.touches(merged_lines_gdf.geometry[connecting_idx])].index)
            for touching_idx in immediately_touching:
                if touching_idx in real_highway_lineids:
                    touches_a_line_thats_already_included = True
                if not touching_idx in contiguous_ids:
                    contiguous_termini += list(merged_lines_gdf.loc[touching_idx,'geometry'].boundary.geoms)
                    contiguous_ids.append(touching_idx)
        if touches_a_line_thats_already_included:
            real_highway_lineids.update(contiguous_ids) 
            continue
        
        long_enough = False
        for terminus in contiguous_termini:
            if not long_enough:
                if terminus.distance(this_segment_termini[0]) > min_length or terminus.distance(this_segment_termini[1]) > min_length:
                    long_enough = True
        if long_enough:
            real_highway_lineids.update(contiguous_ids) 
            continue
        
    real_highways_gdf_utm = merged_lines_gdf.loc[list(real_highway_lineids)]
    real_highways_gdf_ll = real_highways_gdf_utm.to_crs(4326)
    end = timeit.default_timer()
    print('Time prep final4: ', end - start)
    
    del multi_car_G
    del merged_lines_gdf
    del separated_major_roads
    gc.collect()
    
    return real_highways_gdf_ll
    
def bbox_from_shp(file_loc):
    with fiona.open(file_loc,'r') as source: 
        bbox = (source[0]['properties']['BBX_LATMN'],
               source[0]['properties']['BBX_LONMN'],
               source[0]['properties']['BBX_LATMX'],
               source[0]['properties']['BBX_LONMX'],)
    #TODO: calculate bounding box from shape instead of referencing variables
    return bbox

queries = {

'healthcare' : '''
[timeout:900][maxsize:1073741824][out:json];
(
node
["amenity"~"hospital|doctors|clinic|pharmacy"]
(poly:"{poly}");
node
["healthcare"~"alternative|birthing_center|centre|midwife|nurse|hospital|doctor|clinic|pharmacy|yes"]
(poly:"{poly}");
way
["amenity"~"hospital|doctors|clinic|pharmacy"]
(poly:"{poly}");
way
["healthcare"~"alternative|birthing_center|centre|midwife|nurse|hospital|doctor|clinic|pharmacy|yes"]
(poly:"{poly}");
);
(._;);
out skel center qt;
''',

'schools' : '''
[timeout:900][maxsize:1073741824][out:json];
(
node
["amenity"~"school|kindergarten"]
(poly:"{poly}");
node
["school"]
 (poly:"{poly}");
way
["amenity"~"school|kindergarten"]
(poly:"{poly}");
way
["school"]
(poly:"{poly}");
rel
["amenity"~"school|kindergarten"]
(poly:"{poly}");
rel
["school"]
(poly:"{poly}");
);
(._;);
out skel center qt;
''',

'libraries' : '''
[timeout:900][maxsize:1073741824][out:json];
(
node
["amenity"~"library|bookcase"]
(poly:"{poly}");
way
["amenity"~"library|bookcase"]
(poly:"{poly}");
);
(._;);
out skel center qt;
''',
}

categories = {
        'healthcare':[],
        'schools':[],
        'libraries':[]
        }

wkbfab = osmium.geom.WKBFactory()
    
class ServiceHandler(osmium.SimpleHandler): #newer
    def __init__(self):
        super().__init__()
        self.locationlist = {
                'healthcare':[],
                'libraries':[],
                'schools':[],
                'bikeshare':[],
                }
        self.carfreelist = []
        

    def node(self, n):
        if 'amenity' in n.tags and n.tags['amenity'] in ['library','bookcase']:
            self.locationlist['libraries'].append((n.location.lon, n.location.lat))
        
        if ( ('amenity' in n.tags and 
               n.tags['amenity'] in ['school','kindergarten']) or
             ('school' in n.tags) ):
            self.locationlist['schools'].append((n.location.lon, n.location.lat))
            
        if ( ('amenity' in n.tags and 
               n.tags['amenity'] in ['hospital','doctors','clinic','pharmacy']) or
             ('healthcare' in n.tags and 
               n.tags['healthcare'] in ['alternative','birthing_center','centre','midwife','nurse','hospital','doctor','clinic','pharmacy','yes']) ):
            self.locationlist['healthcare'].append((n.location.lon, n.location.lat))

        if ( ('amenity' in n.tags and 
               n.tags['amenity'] in ['bicycle_rental']) or
             ('bicycle_rental' in n.tags) ):
            if 'bicycle_rental' in n.tags:
                if not n.tags['bicycle_rental'] in ['shop']:
                    self.locationlist['bikeshare'].append((n.location.lon, n.location.lat))
            else:
                self.locationlist['bikeshare'].append((n.location.lon, n.location.lat))
                

    def area(self, a):
        try:
            if 'amenity' in a.tags and a.tags['amenity'] in ['library','bookcase']:
                    wkb = wkbfab.create_multipolygon(a)
                    poly = shapely.wkb.loads(wkb, hex=True)
                    centroid = poly.representative_point()
                    self.locationlist['libraries'].append((centroid.x, centroid.y))
                
            if 'amenity' in a.tags and a.tags['amenity'] in ['bicycle_rental']:
                    if 'bicycle_rental' in a.tags:
                        if a.tags['bicycle_rental'] in ['shop']:
                            wkb = wkbfab.create_multipolygon(a)
                            poly = shapely.wkb.loads(wkb, hex=True)
                            centroid = poly.representative_point()
                            self.locationlist['bikeshare'].append((centroid.x, centroid.y))
                    else:
                        wkb = wkbfab.create_multipolygon(a)
                        poly = shapely.wkb.loads(wkb, hex=True)
                        centroid = poly.representative_point()
                        self.locationlist['bikeshare'].append((centroid.x, centroid.y))
                
            if ( ('amenity' in a.tags and 
                   a.tags['amenity'] in ['school','kindergarten']) or
                 ('school' in a.tags) ):
                wkb = wkbfab.create_multipolygon(a)
                poly = shapely.wkb.loads(wkb, hex=True)
                centroid = poly.representative_point()
                self.locationlist['schools'].append((centroid.x, centroid.y))
                
            if ( ('amenity' in a.tags and 
                   a.tags['amenity'] in ['hospital','doctors','clinic','pharmacy']) or
                 ('healthcare' in a.tags and 
                   a.tags['healthcare'] in ['alternative','birthing_center','centre','midwife','nurse','hospital','doctor','clinic','pharmacy','yes']) ):
                wkb = wkbfab.create_multipolygon(a)
                poly = shapely.wkb.loads(wkb, hex=True)
                centroid = poly.representative_point()
                self.locationlist['healthcare'].append((centroid.x, centroid.y))
                
            carfree = False
            if 'leisure' in a.tags and a.tags['leisure'] in ['park', 'playground', 'sports_pitch']:
                if 'access' not in a.tags or a.tags['access'] != 'private':
                    carfree = True
        
            # Check for forests and recreation grounds
            if 'landuse' in a.tags and a.tags['landuse'] in ['forest', 'recreation_ground']:
                if 'access' not in a.tags or a.tags['access'] != 'private':
                    carfree = True
            if carfree:
                wkb = wkbfab.create_multipolygon(a)
                poly = shapely.wkb.loads(wkb, hex=True)
                self.carfreelist.append(poly)
                
        except RuntimeError:
            print('RUNTIME ERROR while finding service area')
            
    def way(self, a):
        try:
            carfree = False
            if 'highway' in a.tags and a.tags['highway'] in ['pedestrian', 'path', 'footway', 'cycleway']:
                if 'access' not in a.tags or a.tags['access'] != 'private':
                    if 'footway' not in a.tags or a.tags['footway'] not in ['sidewalk', 'crossing']:
                        carfree = True
            if carfree:
                wkb = wkbfab.create_linestring(a)
                poly = shapely.wkb.loads(wkb, hex=True)
                self.carfreelist.append(poly)
            
        except RuntimeError:
            print('RUNTIME ERROR while finding service way')
        

#haven't used in ages! Need to make sure it returns x,y instead of lat,lon (y,x)
# def get_point_locations(poly, query):
#     #returns a dictionary
    
#     overpass_url = "http://overpass-api.de/api/interpreter" 
#     #check this is good before production
    
#     poly_str = ox.get_polygons_coordinates(poly)[0]
    
#     services = []
    
#     print ('Querying OSM for locations...')
#     data = ox.overpass_request(data={'data':query.format(poly=poly_str)}, timeout=900)
#     for element in data['elements']:
#         if element['type'] == 'node':
#             services.append(
#                     (element['lat'],element['lon'])
#                     )
#         elif 'center' in element:
#             services.append(
#                     (element['center']['lat'],
#                      element['center']['lon'])
#                     )
    
#     return services
    
