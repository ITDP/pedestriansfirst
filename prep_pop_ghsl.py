import pandas as pd
import geopandas as gpd
gpd.options.use_pygeos = False
import osmnx as ox
import numpy as np
import shapely
from shapely.geometry import Polygon
import rasterstats
import cenpy
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
        
def get_tracts_usa(placetype,
                   placenames,
                    bounds_poly_latlon,
                    cenpy_kwargs = {},
                    ):
    acs = cenpy.products.ACS(2019)    
    tract_census_tables = [
        'B01003',# Total Population
        #'B08201',# Household Size by Vehicles Available
            #TODO: use this to do a more sophisticated estimate of vehicle availability
        'B25046',# Aggregate Number of Vehicles Available by Tenure
        'B19081',# Mean Household Income of Quintiles
        ]
    
    all_tracts = []
    for placename in placenames:
        if placetype == 'place':
            all_tracts.append(acs.from_place(placename, level='tract', variables=tract_census_tables, strict_within=False))
        if placetype == 'msa':
            all_tracts.append(acs.from_msa(placename, level='tract', variables=tract_census_tables, strict_within=False))
    tracts = pd.concat(all_tracts, ignore_index=True)
    tracts_ll = tracts.to_crs(4326)
    
    tracts_clipped = tracts_ll.intersection(bounds_poly_latlon)
    selection = ~tracts_clipped.geometry.is_empty
    select_tracts = tracts_ll[selection]
    select_tracts_utm = ox.project_gdf(select_tracts)
    
    select_tracts_utm['population'] = select_tracts_utm['B01003_001E']
    select_tracts_utm['perc_w_cars'] = select_tracts_utm['B25046_001E'] / select_tracts_utm['population']
    
    for quint in range(1,6):
        select_tracts_utm[f'quint_{quint}_income'] = select_tracts_utm[f'B19081_00{quint}E']
    select_tracts_utm['population_per_sqm'] = select_tracts_utm['population'] / select_tracts_utm.area
    
    # drop unncessary columns (census cols)
    for column in select_tracts_utm.columns:
        if column[:6] in tract_census_tables:
            select_tracts_utm.drop([column], axis=1) #TODO MAKE INPLACE
    select_tracts_utm.drop(['NAME',
                            'state',
                            'county',
                            'tract'], axis=1)
    
    select_tracts = select_tracts_utm.to_crs(4326)
    select_tracts['id'] = select_tracts['GEOID']
    for bike_lts in cyclist_distribution.keys():
        select_tracts[bike_lts] = cyclist_distribution[bike_lts]
    
    return select_tracts
    
    
    #TODO - use maup to allow for conversion to arbitrary geometry (eg grid)
    #Requires higher-resolution source of population data than tracts - eg GHSL
    #Because tracts and grid are both on about the same order of size, so all you do is lose resolution
    # tract_pieces = maup.intersections(tracts, grid, area_cutoff=0)
    # import pdb; pdb.set_trace()
    # tract_pieces_weights = tracts['population']#.groupby(maup.assign(tracts_pieces, grid)).sum()
    # tract_pieces_weights = maup.normalize(tract_pieces_weights, level=0)
    # for variable in ['population',
    #                  'perc_w_cars',
    #                  'quint_1_income',
    #                  'quint_2_income',
    #                  'quint_3_income',
    #                  'quint_4_income',
    #                  'quint_5_income',
    #                  ]:
    #     grid[variable] = maup.prorate(
    #         tract_pieces,
    #         tracts[variable],
    #         weights=tract_pieces_weights,
    #         )
    

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

