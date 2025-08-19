import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = False
import osmnx as ox
import numpy as np
import shapely
from shapely.geometry import Polygon
import rasterstats
import maup
#import utm

# see http://www1.coe.neu.edu/~pfurth/Other%20papers/Dill%202013%204%20types%20of%20cyclists%20TRR.pdf
cyclist_distribution = {
    'bike_lts4': 0.04,
    'bike_lts2': 0.09,
    'bike_lts1': 0.56,
    'bike_lts0': 0.31, #will not bike
    }

def build_grid(city_crs, bounds_poly_utm, low_resolution, exception_gdf_utm=None, high_resolution=None):
    xmin,ymin,xmax,ymax =  bounds_poly_utm.bounds
    # thank you Faraz (https://gis.stackexchange.com/questions/269243/creating-polygon-grid-using-geopandas)
    rows = int(np.ceil((ymax-ymin) /  low_resolution))
    cols = int(np.ceil((xmax-xmin) / low_resolution))
    XleftOrigin = xmin
    XrightOrigin = xmin + low_resolution
    YtopOrigin = ymax
    YbottomOrigin = ymax - low_resolution
    lowres_cells = []
    exception_cells = []
    for i in range(cols):
        Ytop = YtopOrigin
        Ybottom = YbottomOrigin
        for j in range(rows):
            cell = Polygon([(XleftOrigin, Ytop), (XrightOrigin, Ytop), (XrightOrigin, Ybottom), (XleftOrigin, Ybottom)])
            cell = cell.intersection(bounds_poly_utm)
            if exception_gdf_utm is not None:
                if not cell.intersects(exception_gdf_utm.unary_union):
                    lowres_cells.append(cell)
                else:
                    exception_cells.append(cell)
            else:
                lowres_cells.append(cell)
            Ytop = Ytop - low_resolution
            Ybottom = Ybottom - low_resolution
        XleftOrigin = XleftOrigin + low_resolution
        XrightOrigin = XrightOrigin + low_resolution
    highres_cells = []
    if exception_gdf_utm is not None:
        for exception_cell in exception_cells:
            highres_cells += build_grid(exception_cell, high_resolution)
    grid = gpd.GeoDataFrame(geometry=lowres_cells + highres_cells, crs=city_crs)
    grid = grid[grid.area > 0]
    grid.reset_index(inplace=True)
    grid.drop('index', axis=1, inplace=True)
    return grid

def build_grid_from_latlon(bounds_or_poly, low_resolution, exception_gdf_utm=None, high_resolution=None):
    if type(bounds_or_poly) == type(['l','i','s','t']):
        bounds_ll = shapely.geometry.box(*bounds_or_poly)
        df = gpd.GeoDataFrame(geometry=[bounds_ll], crs=4326)
    else:
        df = gpd.GeoDataFrame(geometry=[bounds_or_poly], crs=4326)
    df_utm = ox.project_gdf(df)
    city_crs = df_utm.crs
    return build_grid(city_crs, df_utm.unary_union, low_resolution, exception_gdf_utm, high_resolution)
        
    

def populate_grid_ghsl(grid, 
                  ghsl_fileloc,
                  ghsl_crs,
                  car_distribution = { #default values are estimates for Nairobi
                      'pct_carfree' : 0.8,
                      'pct_onecar': 0.15,
                      'pct_twopluscars': 0.05,
                      },
                  ):
    grid['area_m2'] = grid.geometry.area
    grid_gdf_proj= grid.to_crs(ghsl_crs)
    pop_densities = rasterstats.zonal_stats(grid_gdf_proj, 
                                            ghsl_fileloc, 
                                            stats=['mean'])
    for idx in grid_gdf_proj.index:
        if pop_densities[idx]['mean'] == None:
            mean_dens = 0
        else:
            mean_dens = pop_densities[idx]['mean']
        mean_dens_per_m2 = (mean_dens / 0.0625) / 1000000
        grid_gdf_proj.loc[idx,'pop_dens'] = mean_dens_per_m2
        grid_gdf_proj.loc[idx,'population'] = mean_dens_per_m2 * grid_gdf_proj.loc[idx,'area_m2']
    
    for bike_lts in cyclist_distribution.keys():
        grid_gdf_proj[bike_lts] = cyclist_distribution[bike_lts]
    
    #TODO: remove cells with 0 population (eg water)
    grid_gdf_latlon = grid_gdf_proj.to_crs(4326)
    return grid_gdf_latlon


#for some reason, the total pop of the grid may be different from the total pop of the area
#I must have made some mistake above
#but for now, this will even things out
#TODO figure out what causes this discrepancy
def adjust_population(grid_gdf_latlon, ghsl_crs, ghsl_fileloc):
    grid_gdf_proj = grid_gdf_latlon.to_crs(ghsl_crs)
    pop_sum = rasterstats.zonal_stats([grid_gdf_proj.unary_union], 
                                            ghsl_fileloc, 
                                            stats='sum')
    total_pop=pop_sum[0]['sum']
    df_sum = grid_gdf_latlon.population.sum()
    print(total_pop, df_sum)
    ratio = total_pop / df_sum
    grid_gdf_latlon.population = grid_gdf_latlon.population * ratio
    return grid_gdf_latlon    
    
def setup_grid_ghsl(bounds_poly_latlon, 
               low_resolution, 
               ghsl_fileloc,
               ghsl_crs,
               adjust_pop = True, #False if you don't want to adjust the populations
               car_distribution = { #default values are estimates for Brazil
                   'pct_carfree' : 0.8,
                   'pct_onecar': 0.2,
                   'pct_twopluscars': 0,
                   },
               ):
    bounds_gdf_latlon = gpd.GeoDataFrame(geometry=[bounds_poly_latlon],crs=4326)
    bounds_gdf_utm = ox.project_gdf(bounds_gdf_latlon)
    city_crs = bounds_gdf_utm.crs
    grid_utm = build_grid(city_crs, bounds_gdf_utm.unary_union, low_resolution)
    grid_gdf_latlon = populate_grid_ghsl(grid_utm, ghsl_fileloc, ghsl_crs)
    if adjust_pop:
        grid_gdf_latlon = adjust_population(grid_gdf_latlon, ghsl_crs, ghsl_fileloc)
    return grid_gdf_latlon

#TODO remove? If I'm just going to be calling from setup.py anyway.
def grid_from_bbox(
        bbox, #lat-lon, [xmin,ymin,xmax,ymax]
        ghsl_fileloc,
        resolution = 1000,
        out_file_dir = 'scenarios/existing_conditions/',
        car_distribution = { #default values are estimates for Brazil
            'pct_carfree' : 0.8,
            'pct_onecar': 0.2,
            'pct_twopluscars': 0,
            },
        ):
    bounds_ll = shapely.geometry.box(*bbox)
    df = gpd.GeoDataFrame(geometry=[bounds_ll], crs=4326)
    df_utm = ox.project_gdf(df)
    grid_ll = setup_grid_ghsl(df_utm.crs,
               df_utm.geometry.unary_union,
               resolution,
               ghsl_fileloc,
               adjust_pop = True,
               car_distribution = car_distribution,
               )
    return grid_ll

