'''Functions to calculate isochrones along networks

Using code developed by ITDP-Brazil, available at
https://colab.research.google.com/drive/1xpzFkH7MPwmtFiIfBP2Lg2yqtvzhiF23
'''



import networkx as nx # Pacote de redes do Python
import osmnx as ox # Pacote OSM Network X
import geopandas as gpd # Pacote Pandas Georreferenciado
import shapely
from shapely.geometry import Point, LineString, MultiLineString, Polygon, MultiPolygon, mapping # Pacote Shapely (apenas partes do pacote)
from shapely.ops import unary_union

from copy import deepcopy # duplicar memória do Python
from tqdm import tqdm # Pacote para visualizar downloads
import pandas as pd 

import pdb


def download_graph(coordinates, distances):
    """
    Criação do grafo de ruas do OSM a partir das coordenadas solicitadas
    """
    max_distance = max(distances)
    
    G = False
    print('Fetching street network')
    for coordinate in tqdm(coordinates, desc='Downloading'):

        if G: # "soma" (merge) com grafo já existente (deepcopy utilizado para não perder grafo entre iterações)
            G = nx.compose(deepcopy(G), ox.graph_from_point(coordinate, 
                                                  distance=max_distance+100,
                                                  network_type='walk'))
        else: # inicializa grafo a partir de todos pontos
            G = ox.graph_from_point(coordinate, distance=max_distance+100, 
                                    network_type='walk')            
    
    return G

def network_buffer(G, distance=100):
    all_gs = gpd.GeoDataFrame(ox.graph_to_gdfs(G,nodes=False,node_geometry=False,fill_edge_geometry=True)).buffer(100)
    new_iso = gpd.GeoSeries(all_gs).unary_union # junção de pontos e edges, mas com buracos nas quadras
    return new_iso

def cut(line, distance):
    # Cuts a line in two at a distance from its starting point, 
    # returns the line up to that distance
    if distance <= 0.0 or distance >= line.length:
        return 
    coords = list(line.coords)
    for i, p in enumerate(coords):
        dist_to_point = line.project(Point(p))
        if dist_to_point == distance:
            return LineString(coords[:i+1])
        if dist_to_point > distance:
            cp = line.interpolate(distance)
            return LineString(coords[:i] + [(cp.x, cp.y)])

def proper_iso_polys(G, center_nodes1, distance=500,
                     buffer=100, infill=2500): #infill is the max size that will be infilled
    """
    Returns Shapely geometry representing network buffer from a list of nodes.
    Includes buffering part-way along edges.
    """
    failures = 0
    polygons = []
    complete_nodepairs = set()
    complete_edge_linestrings = []
    partial_edge_linestrings = []
    
    for center_node in tqdm(center_nodes1):
        try:
            distances, routes = nx.single_source_dijkstra(
                G, 
                center_node, 
                target = None, 
                cutoff = distance, 
                weight = 'length')
            # intermediate_nodes = set()
            intermediate_nodes = list(distances.keys())
            for route in routes.values():
                # intermediate_nodes.update(route[1:-1])
                complete_nodepairs.update([(route[i],route[i+1]) for i in range(0,len(route)-1)])
            for origin_node in routes.keys():
                remaining_distance = distance - distances[origin_node]
                for target_node in G[origin_node]:
                    if target_node not in intermediate_nodes:
                        for edge_id in G[origin_node][target_node]:
                            edge = G[origin_node][target_node][edge_id]
                            if 'geometry' in edge:
                                linestring = edge['geometry']
                            else:
                                n_fr = G.nodes[origin_node]
                                n_to = G.nodes[target_node]
                                linestring = LineString([(n_fr['x'],n_fr['y']),(n_to['x'],n_to['y'])])
                            cut_line = cut(linestring, remaining_distance)
                            if cut_line is not None:
                                partial_edge_linestrings.append(cut_line)
        except nx.exception.NodeNotFound:
            pass
    for nodepair in complete_nodepairs:
        for edge_id in G[nodepair[0]][nodepair[1]]:
            edge = G[nodepair[0]][nodepair[1]][edge_id]
            if 'geometry' in edge:
                linestring = edge['geometry']
            else:
                n_fr = G.nodes[nodepair[0]]
                n_to = G.nodes[nodepair[1]]
                linestring = LineString([(n_fr['x'],n_fr['y']),(n_to['x'],n_to['y'])])
            complete_edge_linestrings.append(linestring)
    
    # Combine complete and partial edge geometries into a single GeoSeries
    all_linestrings = gpd.GeoSeries(complete_edge_linestrings + partial_edge_linestrings)
    
    # Apply a buffer around all linestrings and combine into a single geometry
    total_area = all_linestrings.buffer(buffer).unary_union 
    
    # Separate the resulting geometry into individual polygons
    if type(total_area) == shapely.geometry.MultiPolygon:
        polys = list(total_area.geoms)
    elif type(total_area) == shapely.geometry.Polygon:
        polys = [total_area]
    else:
        polys = []
        
        
    # Fill holes in polygons if their area is below the specified threshold    
    infills = []
    for poly in polys:
        for hole in poly.interiors:
            hole_poly = Polygon(hole)
            if hole_poly.area <= infill:
                infills.append(hole_poly)
                
    # Combine the polygons and filled holes into the final geometry
    total_area = gpd.GeoSeries(polys+infills).unary_union
    
    # # output a gdf instead
    # gdf = gpd.GeoDataFrame({
    # "id": list(center_nodes1),
    # "geometry": list(polys+infills)
    # })
    
    
    
    return total_area
    

# G = G_allroads
# distance=500
# edge_buff=100
# node_buff=0
# infill=True
# union = True

# center_nodes1 = center_nodes["pnab"]


def make_iso_polys(G, center_nodes1, distance=500, edge_buff=100, node_buff=0, infill=True, union = True):
    """
    Returns Shapely geometry representing network buffer from a list of nodes.
    """
    failures = 0
    polygons = []
    
    for center_node in tqdm(center_nodes1):
        subgraph = nx.ego_graph(G, center_node, radius=distance, 
                                distance='length')
        
        node_points = [Point((data['x'], data['y'])) 
                       for node, data in subgraph.nodes(data=True)]
        #I used to try and except and keyerror here
        nodes_gdf = gpd.GeoDataFrame({'id': subgraph.nodes()},
                                     geometry=node_points)
    
        nodes_gdf = nodes_gdf.set_index('id')

        edge_lines = []
        for n_fr, n_to in subgraph.edges():
            f = nodes_gdf.loc[n_fr].geometry
            t = nodes_gdf.loc[n_to].geometry
            edge_lines.append(LineString([f,t]))

        n = nodes_gdf.buffer(node_buff).geometry
        e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry
        all_gs = list(n) + list(e)
        new_iso = gpd.GeoSeries(all_gs).unary_union 
        
        if infill:
            try:# try to fill in surrounded areas so shapes will appear solid
                new_iso = Polygon(new_iso.exterior)
            except AttributeError: #empty geometrycollection?
                pass
        polygons.append(new_iso)
        
    # isochrone_polys = unary_union(polygons)
      
    # output a gdf
    gdf = gpd.GeoDataFrame({
    "id": list(center_nodes1),
    "geometry": list(polygons)
    })
    
    # junta todos os poligonos no entorno das coordenadas 

    return gdf, failures



