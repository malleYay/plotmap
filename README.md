# plotmap
generate large-scale maps for pen plotters from OSM/shapefile/geotiff data

## Installation
- install python-packages.txt
- download data
  - bathymetric data in GeoTIFF-format from https://www.gebco.net/data_and_products/gridded_bathymetry_data/
  - coastlines (WGS84 projection) in shapefile-format from https://osmdata.openstreetmap.de/data/coastlines.html
  - land and water polygons as shapefile (SHP) from https://www.naturalearthdata.com/downloads/
- copy data in world_data
- adjust setup block in new_map_world.py
- run new_map_world.py
