import subprocess
import fiona
import os
import json

import people_near_services


def from_id_hdc(hdc, folder = None):
    #select city from ID number
    with fiona.open('GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp','r') as ucdb:
        for city in ucdb:
            if city['properties']['ID_HDC_G0'] == hdc:
                test_city = city
    #save city geometry so that I can take an extract from planet.pbf within it
    if not os.path.isdir(str(hdc)):
        os.mkdir(str(hdc))
    with open(str(hdc)+'/boundaries.geojson', 'w') as out:
        out.write(json.dumps(test_city))
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
    
    if not folder:
        folder = str(hdc)
    results = people_near_services.pnservices(test_city, folder_name = folder)
    print(str(results))

hdcs = {
'Beijing': 10687,
'Tianjin': 10922,
'Delhi': 6955,
'Karachi': 6169,
'Guangzhou': 12080,
'Mexico City': 154,
'Pune': 7041,
'Manila': 12829,
'Chennai': 8675,
'Bogota': 621,
'Medan': 10692,
'Nairobi': 4808,
'Jakarta': 11862,
'Brasilia': 1210,
'Rio de Janeiro': 1361,
'Buenos Aires': 1105,
'Dar es Salaam': 5222,
'Ahmadabad': 6651
        }

hdcs = {'Havana':473}

for city in hdcs.keys():
    from_id_hdc(hdcs[city])
