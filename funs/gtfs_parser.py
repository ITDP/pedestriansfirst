import requests
import csv
import os

import gtfs_kit as gk
import zipfile
import wget
import io

import subprocess
import datetime
import shutil
import string
import os

import pdb
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Point
import re
import math

from zipfile import ZipFile

from funs.get_jurisdictions import get_jurisdictions

# filename = "input_data/gtfs/2024/gtfs_08154_202412_2010Etufor.zip"
# hdc = "08154"

def feed_from_filename(filename, hdc):
    
    try:
        if not os.path.isdir(f'temp_gtfs_dir/{hdc}/'):
            # os.mkdir('temp_gtfs_dir')
            os.mkdir(f'temp_gtfs_dir/{hdc}/')
        with ZipFile(filename, 'r') as zgtfs:
            zgtfs.extractall(f'temp_gtfs_dir/{hdc}/')
        # command = 'unzip '+filename+' -d temp_gtfs_dir/'
        # print(command)
        # subprocess.check_call(command, shell=True)
    except:
        if os.path.exists(f'temp_gtfs_dir/{hdc}/'):
            shutil.rmtree(f'temp_gtfs_dir/{hdc}/')
        return False
    if os.path.exists(f'temp_gtfs_dir/{hdc}/calendar.txt'):
        #this fixes a bug that was happening because San Francisco, of all places, ended text lines with whitespace
        with open(f'temp_gtfs_dir/{hdc}/calendar.txt','r') as calfile:
            out = ''
            for line in calfile:
                out += line.strip()
                out += '\n'
        with open(f'temp_gtfs_dir/{hdc}/calendar.txt','w') as calfile:
            calfile.write(out)
    try:
        feed = gk.read_feed(f'temp_gtfs_dir/{hdc}/', dist_units = 'km')
    except pd.errors.ParserError:
        return False
    except KeyError:
        return False
    if os.path.exists(f'temp_gtfs_dir/{hdc}/'):
        shutil.rmtree(f'temp_gtfs_dir/{hdc}/')
    return feed

def log(folder_name, msg):
    print(msg)
    with open(folder_name+"/temp/gtfs/log.txt", "a") as myfile:
        myfile.write(msg)


def get_stop_frequencies(feed, headwaylim, folder_name, filename):
    #if "error" in validation['type']: #turns out this removes files that would otherwise succeed
        #remove the gtfs file so it doesn't cause problems for r5py
   #     log(folder_name, "validation_failed,"+feed.agency.agency_name[0])
    #    os.remove(filename)
     #   return {}
     # feed = feed_from_filename("input_data/gtfs/2023/gtfs_05472_202404_Transjakarta_fixed.zip")
     # feed = feed_from_filename("input_data/gtfs/2023/gtfs_05472_202404_Bogor_Angkots__fixed.zip")
     # feed = feed_from_filename("input_data/gtfs/2023/raw/gtfs_05472_202404_Bogor_Angkots_.zip") 
     # feed = feed_from_filename("input_data/gtfs/2023/raw/gtfs_05472_202404_Transjakarta.zip")
     # feed = feed_from_filename("/home/kauebraga/Downloads/gtfs_05472_202404_Transjakarta.zip")
     # feed = feed_from_filename("/home/kauebraga/Downloads/gtfs_05472_202404_Transjakarta (1).zip")
     # feed = feed_from_filename("input_data/gtfs/2023/gtfs_00070_202404_Gruppo_Torinese_Trasporti.zip")
     # feed = feed_from_filename("input_data/gtfs/2023/gtfs_00047_202404_Freetown_PodaPoda_Sierra_Leone_Road_Transport_Corp.zip")
     # feed = feed_from_filename("input_data/gtfs/2023/gtfs_01156_202404_TMT_Trujillo _Gerencia_de_Transporte_fixed.zip")
     
    try:
        days = feed.get_first_week()[0:5]
    except:
        return {}
    counts = gpd.GeoDataFrame(geometry=[], crs=4326)
    try:
        stopstats = gk.stops.compute_stop_stats(feed, days, 
                                      headway_start_time= '07:00:00', 
                                      headway_end_time= '21:00:00', 
                                      split_directions = False)
    except TypeError:
        log(folder_name, "typeerror,"+feed.agency.agency_name[0]+"\n")
        os.remove(filename)
        return {}
    except ValueError:
        log(folder_name, "valueerror,"+feed.agency.agency_name[0]+"\n")
        os.remove(filename)
        return {}
    if stopstats.empty:
        print("did not get counts (stopstats.empty)", feed.agency)
        os.remove(filename)
        return {}
    for stop_id in stopstats.stop_id.unique():
        headway = stopstats.loc[stopstats['stop_id']==stop_id].mean_headway.mean()
        if headway <= headwaylim and headway != 0:
            try:
                row = feed.stops.loc[feed.stops['stop_id']==stop_id].iloc[0]
                lat = row['stop_lat']
                lon = row['stop_lon']
                counts.loc[stop_id,'headway'] = headway
                counts.loc[stop_id,'geometry'] = Point(lon, lat)
            except IndexError:
                log(folder_name, "indexerror,"+feed.agency.agency_name[0]+"\n")
            except ValueError:
                log(folder_name, "valueerror2,"+feed.agency.agency_name[0]+"\n")
    if not counts.empty:
        try:
            log(folder_name,"success,"+feed.agency.agency_name[0]+"\n")
        except AttributeError:
            log(folder_name,"success,"+"Unknown Name"+"\n")
    else:
        log(folder_name,"counts.empty,"+feed.agency.agency_name[0]+"\n")
    return counts
    
def get_frequent_stops(hdc, year, folder_name, headwaylim = 20):
    # filenames = get_GTFS_from_mobility_database(poly, folder_name+'temp/gtfs/')
   #  hdc = "08154"
    # year = 2024
    # headwaylim = 20
    # folder_name = '/media/kauebraga/data/pedestriansfirst/cities_out'+'/ghsl_region_'+str(hdc)+'/'
    
    directory = f"input_data/gtfs/{year}"
    
    
    if not os.path.exists(folder_name+'temp/gtfs/'):
      os.mkdir(folder_name+'temp/gtfs/')

    filenames = [f for f in os.listdir(directory) if re.search(f"gtfs_{hdc}", f)]
    filenames = [os.path.join(directory, f) for f in filenames]
    # filenames = ''.join(str(x) for x in filenames)
    
    all_freq_stops = gpd.GeoDataFrame(geometry=[], crs=4326)
    wednesdays = []
    for filename in filenames:
        # filename = filenames[0]
        try:
            feed = feed_from_filename(filename, hdc)
            counts = get_stop_frequencies(feed, headwaylim, folder_name, filename)
        except UnicodeDecodeError:
            print ('did not add stops!! UnicodeDecodeError')
            return all_freq_stops, wednesdays
        try:
            all_freq_stops = gpd.GeoDataFrame(pd.concat([all_freq_stops, counts], ignore_index=True),crs = 4326)
            wednesdays.append(feed.get_first_week()[2])
        except TypeError:
            print ('did not add stops!! concat typeperror')
    return all_freq_stops, wednesdays
    
