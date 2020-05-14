'''Functions to calculate isochrones along networks

Using code developed by ITDP-Brazil, available at
https://colab.research.google.com/drive/1xpzFkH7MPwmtFiIfBP2Lg2yqtvzhiF23
'''



import networkx as nx # Pacote de redes do Python
import osmnx as ox # Pacote OSM Network X
import geopandas as gpd # Pacote Pandas Georreferenciado
from shapely.geometry import Point, LineString, Polygon, mapping # Pacote Shapely (apenas partes do pacote)
from shapely.ops import cascaded_union
ox.config(log_console=True, use_cache=True)

from copy import deepcopy # duplicar memória do Pyhton
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

def make_iso_polys(G, center_nodes, distance=500,
                   edge_buff=100, node_buff=0, infill=True):
    """
    Returns Shapely geometry representing network buffer from a list of nodes.
    """
    failures = 0
    polygons = []
    
    for i, center_node in enumerate(center_nodes):
        subgraph = nx.ego_graph(G, center_node, radius=distance, 
                                distance='length')
        
        node_points = [Point((data['x'], data['y'])) 
                       for node, data in subgraph.nodes(data=True)]
        
        try:
            nodes_gdf = gpd.GeoDataFrame({'id': subgraph.nodes()},
                                         geometry=node_points)
        
            nodes_gdf = nodes_gdf.set_index('id')

            edge_lines = []
            for n_fr, n_to in subgraph.edges():
                f = nodes_gdf.loc[n_fr].geometry ##.geometry
                t = nodes_gdf.loc[n_to].geometry
                edge_lines.append(LineString([f,t]))

            n = nodes_gdf.buffer(node_buff).geometry
            e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry
            all_gs = list(n) + list(e)
            new_iso = gpd.GeoSeries(all_gs).unary_union 
            
            if infill:
                # try to fill in surrounded areas so shapes will appear solid
                new_iso = Polygon(new_iso.exterior)
            polygons.append(new_iso)
        except KeyError:
            print ("FAILURE AT",i,center_node)
            failures += 1
        
    isochrone_polys = cascaded_union(polygons) 
    # junta todos os poligonos no entorno das coordenadas 

    return isochrone_polys, failures



