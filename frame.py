import subprocess
import fiona
import os
import os.path
import json
import shutil

import pdb

import people_near_services


def from_id_hdc(hdc, folder = None, kwargs = {}):
    #select city from ID number
    with fiona.open('GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp','r') as ucdb:
        for city in ucdb:
            if int(city['properties']['ID_HDC_G0']) == int(hdc):
                target = city
    return from_city(target, kwargs=kwargs)

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
        if os.path.exists('all_results.json'):
            with open('all_results.json','r') as in_file:
                all_results = json.load(in_file)
        else:
            all_results = {}
        if not str(city['properties']['ID_HDC_G0']) in all_results.keys():
            results = from_city(city)
            all_results.update({city['properties']['ID_HDC_G0']:results})
            with open('all_results.json','w') as out_file:
                json.dump(all_results, out_file)

all_cities()

#hdcs = { #test
#'Mexico City': 154,
#        }



#for city in hdcs.keys():
#    if not os.path.exists(str(hdcs[city])+'/results.json'):
#        from_id_hdc(hdcs[city])
#    else:
#        for file in ['city.o5m','cityhighways.o5m','citywalk.o5m']:
#            if os.path.exists(str(hdcs[city])+'/'+file):
#                os.remove(str(hdcs[city])+'/'+file)

