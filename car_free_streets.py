

import osmnx as ox

ped_filter = ('["area"!~"yes"]'
                  '["highway"~"path|steps|pedestrian|footway"]'
                  '["footway"!~"sidewalk|crossing"]'
                  '[!"crossing"][!"sidewalk"]'
                  '["foot"!~"no"]["service"!~"private"]["access"!~"private"]')

park_filter = ('["foot"!~"no"]["service"!~"private"]["access"!~"private"]')
park_infra = 'way["leisure"="park"]'

  

def get_carfree_graph(poly):
    
    ped_G = 0
    park_G = 0
    
    try:
        ped_G = ox.graph_from_polygon(poly, retain_all = True,
                                  custom_filter = ped_filter, simplify=False)
    except ox.EmptyOverpassResponse:
        pass
    except ConnectionError:
        print ("ConnectionError!")
    
    try:
        park_G = ox.graph_from_polygon(poly, retain_all = True, 
                                   infrastructure = park_infra,
                                   custom_filter = park_filter,
                                   simplify = False
                                   )
    except ox.EmptyOverpassResponse:
        pass

    if ped_G and park_G:
        return ped_G, park_G
    elif ped_G:
        return ped_G, 0
    elif park_G:
        return 0, park_G
    else:
        return 0,0
    
    
