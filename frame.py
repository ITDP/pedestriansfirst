import subprocess
import fiona
import os

#import people_near_services


def from_id_hdc(hdc):
    #select city from ID number
    with fiona.open('GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp','r') as ucdb:
        for city in ucdb:
            if city['properties']['ID_HDC_G0'] == hdc:
                test_city = city
        source_crs = ucdb.crs
        source_schema = ucdb.schema
    #save city geometry so that I can take an extract from planet.pbf within it
    os.mkdir(str(hdc))
    with fiona.open(str(hdc)+'/boundaries.geojson',
                    'w',
                    schema = source_schema,
                    crs = source_crs,
                    driver = 'GeoJSON') as out:
        out.write(test_city)
    #take extract from planet.pbf
    command = "osmium extract -p {}/boundaries.geojson -d {}/ -s simple -v -o city.pbf".format(str(hdc), str(hdc))
    subprocess.call(command.split(' '))
    command = "osmconvert {}/city.pbf >{}/city.pbf".format(str(hdc),str(hdc))
    subprocess.call(command.split(' '))
    command = '''osmfilter {}/city.o5m 
    --drop="area=yes highway=link =motor =proposed 
    =construction =abandoned =platform =raceway 
    service=parking_aisle =driveway =private foot=no" 
    --keep="highway" >{}/city.osm
    '''.format(str(hdc),str(hdc))
    subprocess.call(command.split(' '))
    
    
    #results = people_near_services.pnservices(city, folder_name=str(hdc)+'/')
    #print(str(results))
    
from_id_hdc(54)