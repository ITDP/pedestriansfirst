import requests
import csv
import os

import gtfs_kit as gk

import pdb


api_key = "0b073186-cfc0-47e2-959a-a2be054025ff"
overpass_url = "https://api.transitfeeds.com/v1/" 

def get_all_locations():
    query = overpass_url+"getLocations"
    pdb.set_trace()
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
    query = overpass_url+"getLatestFeedVersion
    params = {'key': api_key,
              'feed': feed_id}
    #pdb.set_trace()
    try:
        resp = requests.get(query, params=params, headers={'Accept-Encoding':'identity'})
        with open('temp_gtfs.zip','wb') as temp:
            temp.write(resp.content)
        feed = gk.read_gtfs('temp_gtfs.zip', dist_units = 'km')
    except:
        feed = None
    if os.path.exists('temp_gtfs.zip'):
        os.remove('temp_gtfs.zip')
    return feed

def get_freq_stops(feed, headwaylim = 20):
    #days = feed.get_first_week()[0:5]
    counts = {}
    days = [
            '20191014',
            '20191015',
            '20191016',
            '20191017',
            '20191018']
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
    if stopstats.empty:
        return {}
    for stop_id in stopstats.stop_id.unique():
        headway = stopstats.loc[stopstats['stop_id']==stop_id].mean_headway.mean()
        if headway <= headwaylim and headway != 0:
            row = feed.stops.loc[feed.stops['stop_id']==stop_id].iloc[0]
            lat = row['stop_lat']
            lon = row['stop_lon']
            counts[stop_id] = [headway, lat, lon]
    if counts:
        print ("got counts!")
    return counts
    
def count_all_sources(sources, headwaylim = 20):
    stop_sets = []
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
    return stop_sets
