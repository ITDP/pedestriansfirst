import requests
import json
import fiona
import osmnx as ox

import osmium
import shapely.wkb

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
            if 'leisure' in a.tags and a.tags['leisure'] in ['park', 'playground']:
                if not('foot' in a.tags and a.tags['foot'] == 'no'):
                    if not('service' in a.tags and a.tags['service'] == 'private'):
                            if not('access' in a.tags and a.tags['access'] == 'private'):
                                carfree = True
            if 'highway' in a.tags and a.tags['highway'] == 'pedestrian':
                if not('foot' in a.tags and a.tags['foot'] == 'no'):
                    if not('service' in a.tags and a.tags['service'] == 'private'):
                            if not('access' in a.tags and a.tags['access'] == 'private'):
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
            if 'leisure' in a.tags and a.tags['leisure'] in ['park', 'playground']:
                if not('foot' in a.tags and a.tags['foot'] == 'no'):
                    if not('service' in a.tags and a.tags['service'] == 'private'):
                            if not('access' in a.tags and a.tags['access'] == 'private'):
                                carfree = True
            if 'highway' in a.tags and a.tags['highway'] in ['pedestrian', 'path','steps','footway']:
                if not('foot' in a.tags and a.tags['foot'] == 'no'):
                    if 'crossing' not in a.tags and 'sidewalk' not in a.tags:
                        if not('footway' in a.tags and a.tags['footway'] in ['sidewalk','crossing']):
                            if not('service' in a.tags and a.tags['service'] == 'private'):
                                if not('access' in a.tags and a.tags['access'] == 'private'):
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
    
