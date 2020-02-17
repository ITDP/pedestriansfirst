import subprocess
import fiona
import os
import os.path
import json
import shutil

import people_near_services


def from_id_hdc(hdc, folder = None, kwargs = {}):
    #select city from ID number
    with fiona.open('GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp','r') as ucdb:
        for city in ucdb:
            if city['properties']['ID_HDC_G0'] == hdc:
                city = city
    return from_city(city, kwargs=kwargs)

def from_city(city, kwargs = {}):
    hdc = city['properties']['ID_HDC_G0']
    #save city geometry so that I can take an extract from planet.pbf within it
    if not os.path.isdir(str(hdc)):
        os.mkdir(str(hdc))
    with open(str(hdc)+'/boundaries.geojson', 'w') as out:
        out.write(json.dumps(city))
    #take extract from planet.pbf
    if not os.path.exists('{}/city.pbf'.format(str(hdc))):
        command = "osmium extract planet-latest.osm.pbf -p {}/boundaries.geojson -s complete_ways -v -o {}/city.pbf".format(str(hdc), str(hdc))
        print(command)
        subprocess.check_call(command.split(' '))
    command = "osmconvert {}/city.pbf -o={}/city.o5m".format(str(hdc),str(hdc))
    print(command)
    subprocess.check_call(command.split(' '))
    command = 'osmfilter {}/city.o5m --keep="highway=" -o={}/cityhighways.o5m'.format(str(hdc),str(hdc))
    print(command)
    subprocess.check_call(command, shell=True)
    command = ['osmfilter {}/cityhighways.o5m --drop="area=yes highway=link =motor =proposed =construction =abandoned =platform =raceway service=parking_aisle =driveway =private foot=no" -o={}/citywalk.o5m'.format(str(hdc),str(hdc))]
    print(command)
    subprocess.check_call(command, shell=True)
    
    folder = str(hdc)+'/'
    
    return people_near_services.pnservices(city, folder_name = folder, **kwargs)

def get_pop(city):
    return city['properties']['P15']

#all cities in descending order
    #todo: read from all_results to skip cities i've already done
def all_cities():
    with fiona.open('GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp','r') as ucdb:
        cities = list(ucdb)
    cities.sort(key=get_pop, reverse = True)
    for city in cities:
        results = from_city(city)
        if os.path.exists('all_results.json'):
            with open('all_results.json','r') as in_file:
                all_results = json.load(in_file)
        else:
            all_results = {}
        all_results.update({city['properties']['ID_HDC_G0']:results})
        with open('all_results.json','w') as out_file:
            json.dump(all_results, out_file)

#all_cities()

hdcs = { #test
#'Delhi':6955,
#'Udaipur': 6712,
#'Abu Dhabi': 5909,
#'Pune': 7041,
#'Aqaba':4388,
#'Amman':4408,
#'Cairo':3902,
#'Brasilia': 1210,
#'Dar es Salaam': 5222,
#'Manila': 12829,
#'Chennai': 8675,
#'Bogota': 621,
#'Medan': 10692,
#'Nairobi': 4808,
#'Jakarta': 11862,
#'Rio de Janeiro': 1361,
#'Buenos Aires': 1105,
#'Ahmadabad': 6651,
#'London':1912,
#'Paris':2046,
#'Amsterdam':2167,
#'Rome':2897,
#'Moscow':3675,
#'New York':945,
#'Toronto':875,
#'DC':855,
#'Boston':1022,
#'Atlanta':559,
#'Miami':556,
#'Houston':315,
#'LA':14,
'SF':10,
'Tel Aviv':3409,
'Istanbul':3562,
'Havana':473,
'Beijing': 10687,
'Tianjin': 10922,
'Delhi': 6955,
'Karachi': 6169,
'Guangzhou': 12080,
'Mexico City': 154,
        }



for city in hdcs.keys():
    if not os.path.exists(str(hdcs[city])+'/results.json'):
        print ('gonna get',hdcs[city])
        from_id_hdc(hdcs[city])
    else:
        for file in ['city.o5m','cityhighways.o5m','citywalk.o5m']:
            if os.path.exists(str(hdcs[city])+'/'+file):
                os.remove(str(hdcs[city])+'/'+file)

#for city in hdcs.keys():
#    try:
#        from_id_hdc(hdcs[city])
#    except Exception as e:
#        with open('ERROR'+city+str(hdcs[city])+'.txt', 'w') as errout:
#            errout.write(str(e))
#            errout.close()
