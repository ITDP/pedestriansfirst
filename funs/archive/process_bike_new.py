import warnings
import datetime
import os
import os.path
import gc
import numpy
import math
import statistics
import rasterstats
import rasterio
import rasterio.mask
import subprocess
import json
import traceback
import shutil
from tqdm import tqdm
#import utm_zone
import numpy as np
import osmnx as ox
import networkx as nx
import pandas as pd
import geopandas as gpd
import shapely.geometry
from shapely.geometry import LineString, Point
from shapely.ops import unary_union, transform
import shapely.ops
import topojson

from  funs.isochrones import *
from  funs.get_service_locations import *
from  funs.gtfs_parser import *
from  funs.prep_bike_osm import *
from  funs.prep_pop_ghsl import *
# import access
#import summarize_ttm

import pdb
import time
 



def process_bike(
	edges,
	quilt_allbike,
	quilt_protectedbike,
	center_nodes,
	):

	# Convert graph to data.frame
	# Preprocessing to handle columns
	# edges = ox.graph_to_gdfs(G_allroads, nodes=False)

	edges = edges.drop(columns=[
		'name',
		'width',
		'oneway',
		'reversed',
		'length',
		'lanes',
		'ref',
		'maxspeed',
		'access',
		'bridge',
		'service',
		'junction',
		'area',
		])

	for col in [
		'highway',
		'cycleway',
		'bicycle',
		'cycleway:left',
		'cycleway:right',
		'cycleway:both',
		]:
		if col not in edges.columns:
			edges[col] = ''

	bike_conditions = (edges['highway'] == 'cycleway') \
		| (edges['highway'] == 'path') & (edges['bicycle']
			== 'designated') | edges['cycleway:right'].isin(['track',
			'lane', 'opposite_track']) | edges['cycleway:left'
			].isin(['track', 'lane', 'opposite_track']) \
		| edges['cycleway:both'].isin(['track', 'lane', 'opposite_track'
			]) | edges['cycleway'].isin(['track', 'lane',
			'opposite_track'])

	# Filter only the bike network

	edges = edges[bike_conditions]
	total_allbike = edges

	# Add a column to tag protected bike infrastructure

	total_allbike['is_protected'] = (total_allbike['highway']
			== 'cycleway') | (total_allbike['highway'] == 'path') \
		& (total_allbike['bicycle'] == 'designated') \
		| total_allbike['cycleway:right'].isin(['track',
			'opposite_track']) | total_allbike['cycleway:left'
			].isin(['track', 'opposite_track']) \
		| total_allbike['cycleway:both'].isin(['track', 'opposite_track'
			]) | total_allbike['cycleway'].isin(['track',
			'opposite_track'])

	# Remove isolated segments for bike infrastructure

	total_allbike['in_real_network'] = 'unknown'

	# Removing isolated lanes (all types) in one step

	def identify_real_network(total_allbike, min_radius=1000,
							  max_jump=300):
		total_allbike['in_real_network'] = 'unknown'  # Initialize column if not already present

		for idx in total_allbike.index:
			already_identified = total_allbike.loc[idx,
					'in_real_network'] == 'unknown'

			# Handle case where 'in_real_network' might be a Series

			if isinstance(already_identified, pd.Series):
				already_identified = already_identified.any()

			if already_identified:
				connected_indices = [idx]
				if shapely.minimum_bounding_radius(total_allbike.loc[connected_indices,
						'geometry'].unary_union) > min_radius:
					total_allbike.loc[connected_indices,
							'in_real_network'] = 'yes'
				else:
					for _ in range(10):  # Safety loop to avoid infinite iterations
						connected_network = \
							total_allbike.loc[connected_indices,
								'geometry'].unary_union
						nearby = \
							total_allbike[total_allbike.distance(connected_network)
								< max_jump]

						if 'yes' in nearby['in_real_network'].unique():
							total_allbike.loc[connected_indices,
									'in_real_network'] = 'yes'
							break

						if set(connected_indices) == set(nearby.index):
							if shapely.minimum_bounding_radius(total_allbike.loc[connected_indices,
									'geometry'].unary_union) \
								> min_radius:
								total_allbike.loc[connected_indices,
										'in_real_network'] = 'yes'
							else:
								total_allbike.loc[connected_indices,
										'in_real_network'] = 'no'
							break
						else:
							connected_indices = list(nearby.index)

					# Final check to ensure network meets minimum radius

					if shapely.minimum_bounding_radius(total_allbike.loc[connected_indices,
							'geometry'].unary_union) > min_radius:
						total_allbike.loc[connected_indices,
								'in_real_network'] = 'yes'

		return total_allbike[total_allbike['in_real_network'] == 'yes']

	# Apply function

	total_allbike = identify_real_network(total_allbike)

	# Separate protected and all-bike networks

	total_protectedbike = total_allbike[total_allbike['is_protected']]

	# Adding to quilt and center_nodes

	quilt_allbike = pd.concat([quilt_allbike, total_allbike])
	quilt_protectedbike = pd.concat([quilt_protectedbike,
									total_protectedbike])

	center_nodes['pnpb'] = set(total_protectedbike.index.map(lambda x: \
							   x[0]).tolist()
							   + total_protectedbike.index.map(lambda x: \
							   x[1]).tolist())
	center_nodes['pnab'] = set(total_allbike.index.map(lambda x: \
							   x[0]).tolist()
							   + total_allbike.index.map(lambda x: \
							   x[1]).tolist())

	return (quilt_allbike, quilt_protectedbike, center_nodes)
