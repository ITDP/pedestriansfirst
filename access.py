
#do this before importing r5py
import sys
sys.argv.append(["--max-memory", "90%"])

import r5py
from r5py import TransportNetwork, TravelTimeMatrixComputer
from r5py import TransportMode
from tqdm import tqdm

import datetime
import numpy as np
import os
import pandas as pd

import prep_pop_ghsl
import prep_bike_osm

#for use of r5py in journey gap calculations
def prepare_mode_settings(**kwargs):
    mode_settings = {}
    general_settings = {
        'departure': datetime.datetime(2023,3,15,8,30),
        #'departure_time_window': datetime.timedelta(hours=1), #this is the default
        #'percentiles': [50], #this is the default
        'max_time':datetime.timedelta(hours=2),
        'max_time_walking':datetime.timedelta(hours=2),
        'max_time_cycling':datetime.timedelta(hours=2),
        'max_time_driving':datetime.timedelta(hours=2),
        'speed_walking':4.8,
        'speed_cycling':12.0,
        'max_public_transport_rides':4,
    }
    general_settings.update(kwargs)
    
    walk_settings = general_settings.copy()
    walk_settings.update({
        'transport_modes':[TransportMode.WALK],
        'access_modes':[TransportMode.WALK],
        })
    mode_settings['WALK'] = walk_settings
    
    transit_settings = general_settings.copy()
    transit_settings.update({
        'transport_modes':[TransportMode.TRANSIT],
        'access_modes':[TransportMode.WALK],
         })
    mode_settings['TRANSIT'] = transit_settings
    
    bike_lts1_settings = general_settings.copy()
    bike_lts1_settings.update({
        'transport_modes':[TransportMode.WALK, TransportMode.BICYCLE],
        'access_modes':[TransportMode.WALK, TransportMode.BICYCLE],
        'max_time_walking':datetime.timedelta(minutes=10),
        'speed_walking':4,
        'max_bicycle_traffic_stress':1
        })
    mode_settings['BIKE_LTS1'] = bike_lts1_settings
    
    bike_lts2_settings = general_settings.copy()
    bike_lts2_settings.update({
         'transport_modes':[TransportMode.WALK, TransportMode.BICYCLE],
         'access_modes':[TransportMode.WALK, TransportMode.BICYCLE],
         'max_time_walking':datetime.timedelta(minutes=10),
         'speed_walking':4,
         'max_bicycle_traffic_stress':2
         })
    mode_settings['BIKE_LTS2'] = bike_lts2_settings
    
    bike_lts4_settings = general_settings.copy()
    bike_lts4_settings.update({
        'transport_modes':[TransportMode.WALK, TransportMode.BICYCLE],
        'access_modes':[TransportMode.WALK, TransportMode.BICYCLE],
        'max_time_walking':datetime.timedelta(minutes=10),
        'speed_walking':4,
        'max_bicycle_traffic_stress':4
        })
    mode_settings['BIKE_LTS4'] = bike_lts4_settings
    
    car_settings = general_settings.copy()
    car_settings.update({
        'transport_modes':[TransportMode.CAR],
        'access_modes':[TransportMode.CAR],
        })
    mode_settings['CAR'] = car_settings
    
    return mode_settings

def value_of_cxn(from_pop, to_dests, t_min):
    #see SSTI's Measuring Accessibility, appendix (p.68)
    #rough average of work and non-work
    baseval = from_pop * to_dests
    return baseval * np.e ** (-0.05 * t_min)


def journey_gap_calculations(
            folder_name,
            current_year,
            boundaries_latlon,
            gtfs_filenames,
            gtfs_wednesdays,
            access_resolution = 2000, #m
            min_pop = 2000, #defaults to 500 ppl/km2
            ):
        for file in [folder_name+'temp/access/grid_pop.geojson',
                      folder_name+'temp/access/city_ltstagged.pbf']:
            if os.path.exists(file):
                os.remove(file)
        
        #prep pop -- it would probably be better to do this straight from GHSL, ugh
        grid_gdf_latlon = prep_pop_ghsl.setup_grid_ghsl(
            boundaries_latlon.unary_union, 
            access_resolution, 
            folder_name+"geodata/population/pop_2020.tif", 
            'ESRI:54009', 
            adjust_pop = True
            )
        grid_gdf_latlon['id'] = grid_gdf_latlon.index
        selection = grid_gdf_latlon.population > min_pop
        grid_gdf_latlon = grid_gdf_latlon[selection]
        grid_gdf_latlon.to_file(folder_name+'temp/access/grid_pop.geojson')
        
        #prep osm (add LTS values)
        original_filename = folder_name+"temp/city.pbf"
        biketagged_filename = folder_name+"temp/access/city_ltstagged.pbf"
        prep_bike_osm.add_lts_tags(original_filename, biketagged_filename)
        
        full_gtfs_filenames = [folder_name+'temp/gtfs/'+name for name in gtfs_filenames]
        print(full_gtfs_filenames)
        try:
            transport_network = TransportNetwork(
                biketagged_filename,
                full_gtfs_filenames
                )
        except:
            print('r5py ERROR')
        
        wednesday_mornings = [datetime.datetime.strptime(wed+' 08:30:00', '%Y%m%d %H:%M:%S') for wed in gtfs_wednesdays]
        latest_wednesday = max(wednesday_mornings)
        mode_settings=prepare_mode_settings(departure = latest_wednesday)
        
        points_gdf_latlon = grid_gdf_latlon.copy()
        points_gdf_latlon.geometry = grid_gdf_latlon.centroid
        
        ttms = {}
        try:
            for mode in ['TRANSIT', 'BIKE_LTS1', 'CAR']:#mode_settings.keys():
                print(f'computing ttm for {mode}')
                ttm_computer = TravelTimeMatrixComputer(transport_network, points_gdf_latlon,**mode_settings[mode])
                ttm_long = ttm_computer.compute_travel_times()
                ttm_wide = pd.pivot(ttm_long, index='from_id', columns='to_id', values='travel_time')
                ttms[mode] = ttm_wide
                ttms[mode].to_csv(folder_name+'temp/access/'+mode+'_ttm.csv')
        except:
           print('FAILED to calculate ttm for', mode)
           return False
            
        #3 versions - cumsum, time, value
        print('calculating indicators for journey gaps')
        for origin_id in tqdm(list(grid_gdf_latlon.index)):
            origin_pop = grid_gdf_latlon.loc[origin_id, 'population']
            grid_gdf_latlon.loc[origin_id, 'time_ratio_weighted_sum'] = 0
            if origin_pop > 0:
                #we weight by both the destination pop and by the gravity model
                # factor, by multiplying by both of them
                # see https://docs.google.com/spreadsheets/d/11SpKFZfN-pr3ieftGoyNBPx5lcqY0b0Wg6RBgdevkMY/edit#gid=0 
                weighting_factor_sum = 0
                for dest_id in grid_gdf_latlon.index:
                    dest_pop = grid_gdf_latlon.loc[dest_id, 'population']
                    if dest_pop > 0 and not origin_id == dest_id:
                        car_time = ttms['CAR'].loc[origin_id, dest_id]
                        sustrans_time = min(ttms['TRANSIT'].loc[origin_id, dest_id],ttms['BIKE_LTS1'].loc[origin_id, dest_id])
                        time_ratio = (sustrans_time/car_time) 
                        #we pretend origin is 1 for now, just to get the weighting factor
                        #we'll weight by origin pop when we calculate analysis-area-wide indicators
                        weighting_factor = value_of_cxn(1, dest_pop, sustrans_time) 
                        time_ratio_with_weighting = time_ratio * weighting_factor
                        if not np.isnan(time_ratio_with_weighting):
                            weighting_factor_sum += weighting_factor
                            grid_gdf_latlon.loc[origin_id, 'time_ratio_weighted_sum'] += time_ratio_with_weighting
                journey_gap = grid_gdf_latlon.loc[origin_id, 'time_ratio_weighted_sum'] / weighting_factor_sum
                grid_gdf_latlon.loc[origin_id, 'journey_gap_unweighted'] = journey_gap
                grid_gdf_latlon.loc[origin_id, 'journey_gap_weighted'] = journey_gap * origin_pop
                
        grid_gdf_latlon.to_file(f'{folder_name}geodata/access/grid_pop_evaluated_{current_year}.geojson')
        return True