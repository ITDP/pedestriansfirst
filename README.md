# Pedestrians First

A Python script to calculate indicators of walkability in cities.

This script uses five indicators to measure walkability in any city on Earth. Those indicators are:
  * People Near Services: the percentage of people living within 1000m of both some form of school and some form of healthcare (including pharmacies).
  * People Near Frequent Transit: the percentage of people living within 500m walking distance of frequent transit (6 departures / hour or better).
  * Car-Free Places: the percentage of people living within 100m of some kind of car-free public place, like a park, plaza, or pedestrian street.
  * Block Density: the size of city blocks, representing the directness of potential walking trips
  * Weighted Population Density: the population density experienced by the average resident. In denser cities, services and amenities are more likely to lie within walkable distance.

The results will be hosted in June on our upcoming interactive platform, [pedestriansfirst.itdp.org](pedestriansfirst.itdp.org)

## Requirements

  * Python 3
  * Python geospatial libraries, including
    * [GeoPandas](https://geopandas.org/)
    * [OSMnx](https://github.com/gboeing/osmnx)
    * [Shapely](https://shapely.readthedocs.io/en/stable/manual.html)
    * [Rasterstats](https://pythonhosted.org/rasterstats/)
    * [Rasterio](https://rasterio.readthedocs.io/en/latest/)
    * [gtfs_kit](https://pypi.org/project/gtfs-kit/)
  * OpenStreetMap data processing tools
    * [Osmium tool](https://github.com/osmcode/osmium-tool/blob/master/README.md)
    * [OSMconvert](https://wiki.openstreetmap.org/wiki/Osmconvert)
    * [OSMfilter](https://wiki.openstreetmap.org/wiki/Osmfilter)
  * Data files
    * [the Urban Centre Database](https://ghsl.jrc.ec.europa.eu/ghs_stat_ucdb2015mt_r2019a.php)
    * [the Global Human Settlement Layer population raster](https://drive.google.com/file/d/1CmB6Wl1Id6GOARypzycFIWDJmFbyfRP4/view?usp=sharing) (this is a file that ITDP has prepared using the European Commission data)
    * [planet-latest.osm.pbf](https://planet.openstreetmap.org/)
    
## Project structure

The project follows the structure:

- `0-process_ghsl`: takes the raw downloaded GHSL and create a table that will list all the GHSL available and the countries the OSM data should be downloaded. Should be run every time we have an update on GHSL;
- `1.1-download_osm`: download the OSM country data for the present moment. Should be run once a year when it is required to update the OSM data;
- `1.2-download_gtfs`: download the GTFS city data for the present moment.  Should be run once a year when it is required to update the GTFS data;
- `2-regional_analysis`: calculate the Pedestrian First indicators for a specific year. 
- `funs`: all the functions required to run the regional analysis script.

    
## License

Pedestrians First is released under the MIT license

Development of this project was supported by the [Bernard van Leer Foundation](https://bernardvanleer.org/)

## Contact

For general inquiries, reach out to [mobility@itdp.org](mobility@itdp.org)

For technical cencerns, reach out to [taylor.reich@itdp.org](taylor.reich@itdp.org)
    
