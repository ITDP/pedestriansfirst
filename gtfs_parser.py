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

import pdb
import pandas


api_key = "0b073186-cfc0-47e2-959a-a2be054025ff"
overpass_url = "https://api.transitfeeds.com/v1/" 

def get_all_locations():
    query = overpass_url+"getLocations"
    resp = requests.get(query, params={'key':api_key}, headers={'Accept-Encoding':'identity'})
    data = resp.json()
    
    return data['results']['locations']

def get_relevant_locs(bbox):
    gtfs_locs = get_all_locations()
    out_locs = []
    for i, loc in enumerate(gtfs_locs):
        if bbox[0] < loc['lat'] < bbox[2] and bbox[1] < loc['lng'] < bbox[3]:
            out_locs.append(loc)
    return out_locs

def get_feed_infos(locations):
    feeds = []
    for loc in locations:
        query = overpass_url+"getFeeds"
        params = {'key': api_key,
                  'location':str(loc['id']),
                  'type': 'gtfs'}
        resp = requests.get(query, params=params, headers={'Accept-Encoding':'identity'})
        results = resp.json()['results']
        if 'feeds' in results.keys():
            feeds = feeds + results['feeds']
    return feeds

def feed_from_id(feed_id):
    if os.path.exists('temp_gtfs.zip'):
        os.remove('temp_gtfs.zip')
    query = overpass_url+"getLatestFeedVersion?key="+api_key+"&feed="+feed_id
    print('QUERY',query)
    try:
        wget.download(query, 'temp_gtfs.zip')
        command = 'unzip temp_gtfs.zip -d temp_gtfs_dir/'
        print(command)
        subprocess.check_call(command, shell=True)
    except:
        if os.path.exists('temp_gtfs.zip'):
            os.remove('temp_gtfs.zip')
        if os.path.exists('temp_gtfs_dir'):
            shutil.rmtree('temp_gtfs_dir')
        return False
    if os.path.exists('temp_gtfs_dir/calendar.txt'):
        #this fixes a bug that was happening because San Francisco, of all places, ended text lines with whitespace
        with open('temp_gtfs_dir/calendar.txt','r') as calfile:
            out = ''
            for line in calfile:
                out += line.strip()
                out += '\n'
        with open('temp_gtfs_dir/calendar.txt','w') as calfile:
            calfile.write(out)
    try:
        feed = gk.read_gtfs('temp_gtfs_dir/', dist_units = 'km')
    except pandas.errors.ParserError:
        return False
    if os.path.exists('temp_gtfs.zip'):
        os.remove('temp_gtfs.zip')
    if os.path.exists('temp_gtfs_dir'):
        shutil.rmtree('temp_gtfs_dir')
    return feed

def feed_from_filename(filename):
    try:
        command = 'unzip '+filename+' -d temp_gtfs_dir/'
        print(command)
        subprocess.check_call(command, shell=True)
    except:
        if os.path.exists('temp_gtfs_dir'):
            shutil.rmtree('temp_gtfs_dir')
        return False
    if os.path.exists('temp_gtfs_dir/calendar.txt'):
        #this fixes a bug that was happening because San Francisco, of all places, ended text lines with whitespace
        with open('temp_gtfs_dir/calendar.txt','r') as calfile:
            out = ''
            for line in calfile:
                out += line.strip()
                out += '\n'
        with open('temp_gtfs_dir/calendar.txt','w') as calfile:
            calfile.write(out)
    try:
        feed = gk.read_gtfs('temp_gtfs_dir/', dist_units = 'km')
    except pandas.errors.ParserError:
        return False
    if os.path.exists('temp_gtfs_dir'):
        shutil.rmtree('temp_gtfs_dir')
    return feed

def get_freq_stops(feed, headwaylim = 20):
    try:
        days = feed.get_first_week()[0:5]
    except:
        return {}
    counts = {}
    try:
        stopstats = gk.stops.compute_stop_stats(feed, days, 
                                          headway_start_time= '05:00:00', 
                                          headway_end_time= '21:00:00', 
                                          split_directions = True)
    except ValueError:
        stopstats = gk.stops.compute_stop_stats(feed, days, 
                                          headway_start_time= '05:00:00', 
                                          headway_end_time= '21:00:00', 
                                          split_directions = False)
    except TypeError:
        return {}
    if stopstats.empty:
        return {}
    for stop_id in stopstats.stop_id.unique():
        headway = stopstats.loc[stopstats['stop_id']==stop_id].mean_headway.mean()
        if headway <= headwaylim and headway != 0:
            try:
                row = feed.stops.loc[feed.stops['stop_id']==stop_id].iloc[0]
                lat = row['stop_lat']
                lon = row['stop_lon']
                counts[stop_id] = [headway, lat, lon]
            except IndexError:
                print("YOLO")
    if counts:
        print ("got counts!")
    return counts
    
def count_all_sources(sources, source_type, headwaylim = 20):
    stop_sets = []
    if source_type == "openmobilitydata":
        for source in sources:
            if 'id' in source.keys():
                feed_id = source['id']
                feed = feed_from_id(feed_id)
                if feed:
                    counts = get_freq_stops(feed, headwaylim = headwaylim)
                    if counts:
                        stop_sets.append(counts)
                        print ("succeded (i hope) for", feed_id)
                    else:
                        print ("failed for", feed_id, 'no frequent service')
                else:
                    print ("failed for", feed_id, 'no feed')
    elif source_type == "local_files":
        for source in sources:
            feed = feed_from_filename(source)
            if feed:
                counts = get_freq_stops(feed, headwaylim = headwaylim)
                if counts:
                    stop_sets.append(counts)
                    print ("succeded (i hope) for", source)
                else:
                    print ("failed for", source, 'no frequent service')
            else:
                print ("failed for", source, 'no feed')
        
    return stop_sets
