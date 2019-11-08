import requests
import csv
import os

import gtfs_kit as gk

import pdb


api_key = "0b073186-cfc0-47e2-959a-a2be054025ff"
overpass_url = "https://api.transitfeeds.com/v1/" 

def get_all_locations():
    #returns ?????
    query = overpass_url+"getLocations?key="+api_key
    resp = requests.get(query)
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
        query = overpass_url+"getFeeds?key="+api_key+"&location="+str(loc['id'])+'&type=gtfs'
        results = requests.get(query).json()['results']
        if 'feeds' in results.keys():
            feeds = feeds + results['feeds']
    return feeds

def feed_from_id(feed_id):
    if os.path.exists('temp_gtfs.zip'):
        os.remove('temp_gtfs.zip')
    query = overpass_url+"getLatestFeedVersion?key="+api_key+"&feed="+feed_id
    #pdb.set_trace()
    try:
        r = requests.get(query)
        with open('temp_gtfs.zip','wb') as temp:
            temp.write(r.content)
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

def get_freq_stopsOLD(feed, headwaylim = 20):
    #days = feed.get_first_week()[0:5]
    days = [
            '20191014',
            '20191015',
            '20191016',
            '20191017',
            '20191018']
    sts = gk.stops.compute_stop_stats(feed, days, freq='1H')
    stop_ids = [t[1] for t in list(sts)]
    counts = {}
    for stop_id in stop_ids:
        count = sum([*sts[('num_trips', stop_id)][5:21],
                     *sts[('num_trips', stop_id)][29:45],
                     *sts[('num_trips', stop_id)][53:69],
                     *sts[('num_trips', stop_id)][77:93],
                     *sts[('num_trips', stop_id)][101:117],])
        count = count / 5
        if count > 960 / headwaylim:
            headway = 960/count
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

def count_all_sourcesOLD(sources, headwaylim = 20):
    stop_sets = []
    for source in sources:
        if 'u' in source.keys():
            if 'd' in source['u'].keys():
                url = source['u']['d']
                feed = feed_from_url(url)
                if feed:
                    counts = get_freq_stops(feed_from_url(url), headwaylim = headwaylim)
                    stop_sets.append(counts)
                    print ("succeded (i hope) for", url)
                else:
                    print ("failed for",url)
    return stop_sets
    
def get_all_frequent_stops_OLD(bbox, frequent = 20):
    all_locs = get_all_locations()
    relevant_locs = get_relevant_locs(bbox, all_locs)
    feeds = get_feed_infos(relevant_locs)
    stop_sets = []
    for feed in feeds:
        if 'u' in feed.keys():
            if 'd' in feed['u'].keys():
                url = feed['u']['d']
                batcmd = '''""C:\\Program Files\\R\\R-3.6.1\\bin\\Rscript.exe" "testtidytransit.r" {0} {1}"'''.format(url, frequent)
                print ("batcmd",batcmd)
                result_code = os.system(batcmd + ' > output.txt')
                print(result_code)
                if os.path.exists('output.txt'):
                    fp = open('output.txt', "r")
                    output = fp.read()
                    fp.close()
                    os.remove('output.txt')
                    print(output)
                if result_code == 0 and os.path.exists('service_out.csv'):
                    print('iffin')
                    stops = {}
                    with open ('service_out.csv') as csvfile:
                        stopreader = csv.reader(csvfile)
                        stopreader.__next__() #skip header
                        for row in stopreader:
                            print (row)
                            if (row[1] not in stops.keys()):
                                stops[row[1]] = row
                            elif float(row[6]) < float(stops[row[1]][7]):
                                stops[row[1]] = row
                        stop_sets.append(stops)
                    os.remove('service_out.csv')
    print ("stop_sets")
    print(stop_sets)
    return stop_sets
                    
