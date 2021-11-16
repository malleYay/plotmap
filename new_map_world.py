from datetime import datetime
import os
import pickle
import sys
import csv
import lxml.etree as ET
import shapely
from shapely.wkb import dumps, loads
from shapely import ops
from shapely.prepared import prep
from shapely.geometry import Point, GeometryCollection, MultiLineString, LineString, Polygon, MultiPolygon
from shapely.geometry import shape, box
import fiona
from HersheyFonts import HersheyFonts
from svgwriter import SvgWriter
import maptools

# ----------------------------------------------------------------------------------------------------
# SETUP >

CACHE_DIRECTORY                 = "./world_data/gebco_2021_sub_ice_topo_geotiff/cache"
SEA_LABELS_FILE                 = "./world_data/labels_sea.geojson"
COASTLINE_FILE                  = "./world_data/simplified-land-polygons-complete-3857/simplified_land_polygons.shp"
BORDER_FILE                     = "./world_data/10m_cultural/ne_10m_admin_0_boundary_lines_land.shp"
CITIES_FILE                     = "./world_data/worldcities.geojson"
BATHYMETRY_DIRECTORY            = "./world_data/gebco_2021_sub_ice_topo_geotiff"
BATHYMETRY_FILES                = [
                                    "gebco_2021_sub_ice_topo_n0.0_s-90.0_w-90.0_e0.0.",
                                    "gebco_2021_sub_ice_topo_n0.0_s-90.0_w-180.0_e-90.0.",
                                    "gebco_2021_sub_ice_topo_n0.0_s-90.0_w0.0_e90.0.",
                                    "gebco_2021_sub_ice_topo_n0.0_s-90.0_w90.0_e180.0.",
                                    "gebco_2021_sub_ice_topo_n90.0_s0.0_w-90.0_e0.0.",
                                    "gebco_2021_sub_ice_topo_n90.0_s0.0_w-180.0_e-90.0.",
                                    "gebco_2021_sub_ice_topo_n90.0_s0.0_w0.0_e90.0.",
                                    "gebco_2021_sub_ice_topo_n90.0_s0.0_w90.0_e180.0."
                                ]
TERRAIN_DIRECTORY               = "world_data/gebco_2021_sub_ice_topo_geotiff"
TERRAIN_FILES                   = [
                                    "gebco_2021_sub_ice_topo_n0.0_s-90.0_w-90.0_e0.0.",
                                    "gebco_2021_sub_ice_topo_n0.0_s-90.0_w-180.0_e-90.0.",
                                    "gebco_2021_sub_ice_topo_n0.0_s-90.0_w0.0_e90.0.",
                                    "gebco_2021_sub_ice_topo_n0.0_s-90.0_w90.0_e180.0.",
                                    "gebco_2021_sub_ice_topo_n90.0_s0.0_w-90.0_e0.0.",
                                    "gebco_2021_sub_ice_topo_n90.0_s0.0_w-180.0_e-90.0.",
                                    "gebco_2021_sub_ice_topo_n90.0_s0.0_w0.0_e90.0.",
                                    "gebco_2021_sub_ice_topo_n90.0_s0.0_w90.0_e180.0."
                                ]
# TODO: move sea labels to world_data/
SEA_LABELS                      = [
                                    [[40.5,     -150],  "Nordpazifik"],
                                    [[57,       168],   "Beringmeer"],
                                    [[-40.6,    -148],  "Pazifik"],
                                    # [[24.201523, -93.767498], "Golf von Mexiko"],
                                    [[13,     -78],   "Karibisches Meer"],
                                    [[59.7,     -90],   "Hudson"],
                                    [[58.7,     -88],   "Bay"],
                                    # [[66.004326, -80.944026], "Nordwestpassage"],
                                    [[74,       -71],   "Baffin Bay"],
                                    [[56,       -53],   "Labradorsee"],
                                    [[33,       -42],   "Nordatlantik"],
                                    [[-38,      -30],   "Atlantik"],
                                    [[68,       -15],   "Nordmeer"],
                                    # [[77.532554, -9.455168], "Groenlandsee"],
                                    [[72,       30],    "Barentssee"],
                                    [[74,       63],    "Karasee"],
                                    [[34.2,     13.2],  "Mittelmeer"],
                                    # [[0.754504, -2.021039], "Golf von Guinea"],
                                    [[12,       58.2],  "Arabisches Meer"],
                                    [[-20,      73.5],  "Indischer Ozean"],
                                    # [[11,       82],    "Golf von Bengalen"],
                                    [[-40,      154],   "Tasmanische"],
                                    [[-43,      158],   "See"],
                                    [[-16,      151],   "Korallenmeer"],
                                    [[25.5,     125],   "Chinesisches Meer"],
                                    [[39.3,     130.5], "Jap. Meer"],
                                    [[56,       142],   "Ochotskisches"],
                                    [[53,       146],   "Meer"],
                                    [[73,       160],   "Ostsibirische See"],
                                ]

MAP_CENTER                      = [0, 0]
MAP_SIZE                        = [1500, 1500]                  # unit for data: m / unit for SVG elements: px or mm
VIEWPORT_OFFSET                 = [0, 250]                      # cutout of min and max
TILE_NUMBERS                    = [3, 1]                        # number of tiles [x, y]
VIEWPORT_SIZE                   = [MAP_SIZE[0]-2*VIEWPORT_OFFSET[0], MAP_SIZE[1]-2*VIEWPORT_OFFSET[1]]
TILE_SIZE                       = [VIEWPORT_SIZE[0]/TILE_NUMBERS[0], VIEWPORT_SIZE[1]/TILE_NUMBERS[1]]
MAP_FRAGMENT_OFFSET             = VIEWPORT_OFFSET
MAP_FRAGMENT_SIZE               = VIEWPORT_SIZE
HOLE_DIST                       = 15                            # distance of mounting holes to the edges

# instead of creating svg(s) of the whole map, create a small svg-preview
PREVIEW                          = True
PREVIEW_SIZE                     = [100, 100]                        # preview size in mm
PREVIEW_POSITION                 = [750, 250]                        # position on map (considering MAP_SIZE and VIEWPORT_OFFSET) in mm

THRESHOLD_CITY_POPULATION       = 1

CONNECT_DATABASE                = False

FONT_SIZE                       = 5
FONT_SIZE_LARGE                 = 12

CREATE_CACHE                    = False

DRAW_META                       = True
DRAW_SCREWHOLES                 = True
DRAW_LAT_LON_LINES              = True
DRAW_SEA_LABELS                 = True
DRAW_COASTLINES                 = True
DRAW_PLACES                     = False
DRAW_CITIES                     = True
DRAW_URBAN_AREAS                = False
DRAW_BORDERS                    = True
DRAW_ADMIN_REGIONS              = False
DRAW_BATHYMETRY                 = True
DRAW_TERRAIN                    = True

CITY_CIRCLE_RADIUS              = 2.0

COLOR_COASTLINES                = [0, 0, 0]
COLOR_TERRAIN                   = [88, 47, 14]
COLOR_BATHYMETRY                = [11, 72, 107]
COLOR_BORDERS                   = [65, 72, 51]
COLOR_SCREWHOLES                = [0, 0, 0]
COLOR_LAT_LON_LINES             = [0, 0, 0]
COLOR_CITIES                    = [0, 0, 0]
COLOR_CITIES_CIRCLES            = [255, 0, 0]
COLOR_URBAN_AREAS               = [0, 0, 0]
COLOR_ADMIN_REGIONS             = [0, 0, 0]
COLOR_SEA_LABELS                = [0, 0, 0]

BG_COLOR                        = "gray"
DARK_MODE                       = False

# bathymetry
# TODO: remove necessity to indicate num_layer, min_height and max_height
BATHYMETRY_NUM_LAYERS = 9 # 15
BATHYMETRY_MIN_HEIGHT = -9000
BATHYMETRY_MAX_HEIGHT = 0
# TODO: remove necessity to indicate num_layer, min_height and max_height
# terrain
TERRAIN_NUM_LAYERS = 18 # 30
TERRAIN_MIN_HEIGHT = 0
TERRAIN_MAX_HEIGHT = 9000

# < SETUP
# ----------------------------------------------------------------------------------------------------


TIMER_STRING                    = "{:<60s}: {:2.2f}s"

DB_NAME                         = "import"
DB_PREFIX                       = "osm_"

SIMPLIFICATION_MAX_ERROR        = 0.2 #1.0 # 0.2                    # unit in map coordinates (px or mm)

MAP_SIZE_SCALE                  = maptools.EQUATOR/MAP_SIZE[0]      # increase or decrease MAP_SIZE by factor

""" ------------------------------

Projection:

Map is using OSMs Web-Mercator (3857) projection
Natural Earth Shapefiles are encoded in WGS84 (4326)

------------------------------ """


# ----------------------------------------------------------------------------------------------------
# functions

def get_text(font, text):

    lines_raw = font.lines_for_text(text)
    lines_restructured = []
    for (x1, y1), (x2, y2) in lines_raw:
        lines_restructured.append([[x1, y1], [x2, y2]])
    lines = MultiLineString(lines_restructured)

    return lines


def simplify_polygon(polys, min_area=None):

    polys_simplyfied = []
    errors_occured = 0
    skipped = 0

    for i in range(0, len(polys)):
        poly = polys[i]

        # simplified_polygon = simplify_polygon(coastlines[i], epsilon=SIMPLIFICATION_MAX_ERROR)
        poly = poly.simplify(SIMPLIFICATION_MAX_ERROR)

        if not type(poly) is Polygon:
            errors_occured += 1
            continue

        if min_area is not None and poly.area < min_area:
            skipped += 1
            continue

        polys_simplyfied.append(poly)

    return polys_simplyfied, errors_occured, skipped

def simplify_linestring(lines):

    lines_simplyfied = []
    errors_occured = 0
    skipped = 0

    for i in range(0, len(lines)):
        line = lines[i]

        line = line.simplify(SIMPLIFICATION_MAX_ERROR)

        if not type(line) is LineString:
            errors_occured += 1
            continue

        lines_simplyfied.append(line)

    return lines_simplyfied, errors_occured, skipped

def get_poly_layer_from_geojson(filename, min_area=None, latlons_flipped=False):

    layers = {}

    shapefile = fiona.open(filename)

    if latlons_flipped:
        func = conv.convert_wgs_to_map_list_lon_lat
    else:
        func = conv.convert_wgs_to_map_list

    counter = 0
    num_items = len(shapefile)
    for item in shapefile:
        counter += 1
        print("processing item {}/{} in {}".format(counter, num_items, filename), end="\r")
        shapefile_geom = shape(item["geometry"])        
        geom = ops.transform(func, shapefile_geom)
        geom = geom.simplify(SIMPLIFICATION_MAX_ERROR)

        prop = item["properties"]
        layer_number = prop["layer"]
        if layer_number not in layers:
            layers[layer_number] = []

        if type(geom) is Polygon:

            if min_area is not None and geom.area < min_area:
                pass
            else:
                layers[layer_number].append(geom)

        elif type(geom) is MultiPolygon:
            for g in geom.geoms:

                if not geom.is_valid:
                    geom = geom.buffer(0.01)

                if min_area is not None and geom.area < min_area:
                    pass
                else:
                    layers[layer_number].append(g)
        else:
            raise Exception("parsing shapefile: unexpected type: {}".format(geom))

    print("")

    # [
    #   [poly, poly, ...],
    #   [poly, poly, ...],
    # ] 

    # problem: if the lowest layers are not populated, they
    # should not be discard altogether but be empty lists

    flat_list = []
    for i in range(0, max(layers.keys())+1):

        if i not in layers:
            flat_list.append([])
        else:
            flat_list.append(layers[i])

    return flat_list

    # [poly, poly, poly, ...]

    # flat_list = []
    # for sublist in list(layers.values()):
    #     for item in sublist:
    #         flat_list.append(item)

    # return flat_list

def get_polys_from_shapefile(filename, min_area=None, latlons_flipped=False):

    geometries = []
    shapefile = fiona.open(filename)

    if latlons_flipped:
        func = conv.convert_wgs_to_map_list_lon_lat
    else:
        func = conv.convert_wgs_to_map_list

    for item in shapefile:
        shapefile_geom = shape(item['geometry'])
        geom = ops.transform(func, shapefile_geom)
        geom = geom.simplify(SIMPLIFICATION_MAX_ERROR)

        if type(geom) is Polygon:
            if not geom.is_valid:
                geom = geom.buffer(0.01)

            if min_area is not None and geom.area < min_area:
                pass
            else:
                geometries.append(geom)
        elif type(geom) is MultiPolygon:
            for g in geom.geoms:

                if not geom.is_valid:
                    geom = geom.buffer(0.01)

                if min_area is not None and geom.area < min_area:
                    pass
                else:
                    geometries.append(g)
        else:
            raise Exception("parsing shapefile: unexpected type: {}".format(geom))

    return geometries

def get_lines_from_shapefile(filename, latlons_flipped=False):

    geometries = []
    shapefile = fiona.open(filename)

    if latlons_flipped:
        func = conv.convert_wgs_to_map_list_lon_lat
    else:
        func = conv.convert_wgs_to_map_list

    for item in shapefile:
        shapefile_geom = shape(item['geometry'])
        geom = ops.transform(func, shapefile_geom)
        geom = geom.simplify(SIMPLIFICATION_MAX_ERROR)

        if type(geom) is LineString:
            geometries.append(geom)
        elif type(geom) is MultiLineString:
            for g in geom.geoms:
                geometries.append(g)
        else:
            raise Exception("parsing shapefile: unexpected type: {}".format(geom))

    return geometries

def write_geometries_to_file(filename, geoms):

    with open(filename, "wb") as handle:
        pickle.dump(geoms, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_geometries_from_file(filename):

    with open(filename, "rb") as handle:
        return pickle.load(handle)

# returns list of polys
def _smooth_poly(p, factor):

    p = p.buffer(factor).buffer(-factor)

    if type(p) is Polygon:

        if not p.is_valid:
            p = p.buffer(0.01)

        return [p]
    else:
        res = []

        for g in p.geoms:

            if not g.is_valid:
                g = g.buffer(0.01)

            res.append(g)

        return res 


def _polygons_to_linestrings(polygons):
    linestrings = []
    for poly in polygons:

        subpolys = []

        if not poly.is_valid:
            poly = poly.buffer(0.01)
        subpolys.append(poly)

        for p in subpolys:

            outline = p.boundary

            if type(outline) is MultiLineString:
                for g in outline.geoms:
                    linestrings.append(g)
            else:
                linestrings.append(outline)

    return linestrings

def polygons_to_linestrings(polygons, flatten=True):

    linestrings = []

    # polygons is a list of layers
    if len(polygons) > 0 and type(polygons[0]) is list: 
        for layer in polygons:
            if flatten:
                linestrings += _polygons_to_linestrings(layer)
            else:
                linestrings.append(_polygons_to_linestrings(layer))
          
    # polygons is a list of polygons  
    else: 
        linestrings = _polygons_to_linestrings(polygons)

    return linestrings

def _unpack_geometrycollection(m):

    polys = []

    if type(m) is GeometryCollection:
        for g in m.geoms:
            polys += _unpack_geometrycollection(g)
    elif type(m) is MultiPolygon:
        for g in m.geoms:
            polys += _unpack_geometrycollection(g)
    elif type(m) is Polygon:
        polys.append(m)
    else:
        print("unknown geometry: {}".format(type(m)))

    return polys

def polygons_merge_tiles(polys):

    unified_poly = ops.unary_union(polys)
    return _unpack_geometrycollection(unified_poly)

def cut_linestrings(linestrings, cut_polys): # may return [MultiLineString, LineString, ...]

    linestrings_processed = []

    cut_poly = ops.unary_union(cut_polys)
    cut_poly = cut_poly.simplify(SIMPLIFICATION_MAX_ERROR)

    for i in range(0, len(linestrings)):

        if i%20 == 0:
            print("cut linestring {:10}/{:10} ({:5.2f})".format(i, len(linestrings), (i/len(linestrings))*100) , end="\r")
        
        line = linestrings[i]
        line = line.simplify(SIMPLIFICATION_MAX_ERROR)

        for l in validate_linestring(line):   
            linestrings_processed.append(l.difference(cut_poly))

    return linestrings_processed

def cut_polygons(polygons, cut_polys):

    polygons_processed = []

    for i in range(0, len(polygons.geoms)):

        if i%20 == 0:
            print("cut polygon {:10}/{:10} ({:5.2f})".format(i, len(polygons.geoms), (i/len(polygons.geoms))*100), end="\r")
        
        poly = polygons.geoms[i]
        poly = poly.simplify(SIMPLIFICATION_MAX_ERROR)
        poly_bounds = poly.bounds
        
        for cut_poly in cut_polys:    

            if maptools.check_polygons_for_overlap(poly_bounds, cut_poly.bounds):
                poly = poly.difference(cut_poly)
        
        polygons_processed.append(poly)

    print("") # newline to counter \r

    return polygons_processed

def recalculate_exclusion_zones(zones):

    zones_simplified = []

    for z in zones:
        zones_simplified.append(z.simplify(SIMPLIFICATION_MAX_ERROR))

    zones_simplified = [ops.unary_union(zones_simplified)]

    # TODO: check for validity

    return zones_simplified

def validate_polygon(poly):

    p = poly

    if p.area <= 1:
        return []

    if not p.is_valid:

        p = p.buffer(0.01)

        if not p.is_valid:
            print("error: polygon not valid")
            return []

    if type(p) is MultiPolygon:
        return list(p.geoms)
    else:
        return [p]

def validate_linestring(line):

    l = line

    if not l.is_valid:
        # print("error: LineString not valid")
        return []

    if l.is_empty:
        # print("error: LineString is empty")
        return []

    if not len(l.bounds) == 4:
        # print("error: LineString is empty (bounds)")
        return []

    if type(l) is MultiLineString:
        lines = []
        for g in l.geoms:
            lines += validate_linestring(g)
        return lines
    elif type(l) is GeometryCollection:

        for g in l.geoms:
            lines = []
            for h in l.geoms:
                lines += validate_linestring(h)
            return lines
    else:
        return [l]

def add_layers_to_writer(poly_layers, svg_handler, cut_polys=[]):

    options = {
        "stroke_width": 0, #0.5
        "opacity": 0,
        "layer": "bathymetry",
        "stroke": [11, 72, 107] # [0, 198, 189]
    }

    cut_poly = ops.unary_union(cut_polys)

    if DARK_MODE:
        options["stroke"] = [255, 255, 255]

    for layer_i in range(0, len(poly_layers)):

        polys = poly_layers[layer_i]
        hatching_name = "bathymetry_hatching_{}".format(layer_i)

        # if layer_i > 8:
        #     options["stroke"] = [11, 72, 107]

        for i in range(0, len(polys)):

            print("\33[2K   adding layer {}/{} | poly {}/{}".format(layer_i, len(poly_layers), i, len(polys)), end="\r")

            p = polys[i]
            p = p.difference(cut_poly)

            save_polys = validate_polygon(p)

            for save_poly in save_polys:
                
                save_poly = save_poly.intersection(viewport_polygon)
      
                for poly in validate_polygon(save_poly):
                    svg_handler.add_polygon(poly, **options, hatching=hatching_name)
                    # svg.add_polygon(poly, **options)

                    # dist = 1 + 0.2*layer_i

                    # outline_polys = [poly]
                    # for d in range(0, 100):
                    #     for p in outline_polys: 
                    #         svg.add_polygon(p, **options)
                        
                    #     new_polys = []
                    #     for p in outline_polys:
                    #         p_smaller = p.buffer(-dist).simplify(SIMPLIFICATION_MAX_ERROR)
                    #         p_smaller_validated = validate_polygon(p_smaller)

                    #         for pv in p_smaller_validated:
                    #             if pv.area > 1.0:
                    #                 new_polys.append(pv)

                    #     outline_polys = new_polys
                    #     if len(outline_polys) == 0:
                    #         break

def cut_bathymetry_inplace(layers, tool_polys):

    current_poly = None
    tool_poly = None

    for i in range(0, len(layers)):
        layer = layers[i]
        for j in range(0, len(layer)):
            current_poly = layer[j]

            if current_poly.area < 1.0:
                continue

            print("cut bathymetry layer {}/{} || geometry {}/{}".format(i, len(layers), j, len(layer)), end="\r")

            for k in range(0, len(tool_polys)):
                tool_poly = tool_polys[k]

                if not current_poly.is_valid:
                    current_poly = current_poly.buffer(0.1)

                if not tool_poly.is_valid:
                    continue

                # if maptools.check_polygons_for_overlap(current_poly, tool_poly):
                current_poly = current_poly.difference(tool_poly)

            layer[j] = current_poly.simplify(SIMPLIFICATION_MAX_ERROR)

    print("")

    return layers

# def cut_bathymetry(layers, tool_polys):

#     for i in range(0, len(layers)):
#         layer = layers[i]

#         mpoly = MultiPolygon(layer)
  
#         for tool_poly in tool_polys:

#             if not mpoly.is_valid:
#                 mpoly = mpoly.buffer(0.1)

#             mpoly = mpoly.difference(tool_poly)

#         layers[i] = mpoly.simplify(SIMPLIFICATION_MAX_ERROR)

#     return layers

def load_bathymetry_file(filename, difference=None):
    cache_file = os.path.join(CACHE_DIRECTORY, filename)

    data = None

    if os.path.exists(cache_file):
        data = load_geometries_from_file(cache_file)
        # if difference is not None:
        #     data = cut_bathymetry_inplace(data, difference)
        print("loaded from cache: {} [{} layers]".format(filename, len(data)))
    else:
        data = get_poly_layer_from_geojson(os.path.join(BATHYMETRY_DIRECTORY, filename), latlons_flipped=True, min_area=5.0)
        if difference is not None:
            data = cut_bathymetry_inplace(data, difference)

        write_geometries_to_file(cache_file, data)
        print("processed and written to cache: {}".format(filename))

    return data

def write_bathymetry_to_cache(filename, difference=None):
    cache_file = os.path.join(CACHE_DIRECTORY, filename)

    data = None

    if not os.path.exists(cache_file):
        data = get_poly_layer_from_geojson(os.path.join(BATHYMETRY_DIRECTORY, filename), latlons_flipped=True, min_area=5.0)
        if difference is not None:
            data = cut_bathymetry_inplace(data, difference)

        write_geometries_to_file(cache_file, data)
        print("processed and written to cache: {}".format(filename))
    else:
        print("skipped writing to cache")

    return

# def read_cities_csv(filename):
#     data = []
#     with open('{}.csv'.format(filename), newline='\n') as csvfile:
#         reader = csv.DictReader(csvfile)
#         for row in reader:
#             data.append([
#                 row['city_ascii'],
#                 row['lat'],
#                 row['lng'],
#                 row['country'],
#                 row['capital'],
#                 row['population']
#                 ])

#     return data

def load_coastline():

    timer_start = datetime.now()

    coastlines = []

    shapefile = fiona.open(COASTLINE_FILE)

    for item in shapefile:
        shp_geom = shape(item["geometry"])

        coastlines.append(shp_geom)

    print(TIMER_STRING.format("loading coastline data", (datetime.now()-timer_start).total_seconds()))   

    timer_start = datetime.now()
    for i in range(0, len(coastlines)):
        coastlines[i] = ops.transform(conv.convert_mercator_to_map_list, coastlines[i])
    print(TIMER_STRING.format("transforming coastline data", (datetime.now()-timer_start).total_seconds()))   

    timer_start = datetime.now()
    coastlines, errors_occured, skipped = simplify_polygon(coastlines, min_area=5.0)

    print(TIMER_STRING.format("simplifing coastline data ({} errors, {} skipped)".format(errors_occured, skipped), (datetime.now()-timer_start).total_seconds())) 

    return coastlines

def draw_meta(svg_handler):

    # options = {
    #     "stroke_width": 2.0,
    #     "layer": "meta",
    #     "stroke": [0, 0, 0],
    #     "opacity": 0
    # }

    options_screwholes = {
        "stroke_width": 2.0,
        "layer": "meta",
        "stroke": COLOR_SCREWHOLES,
        "opacity": 0
    }

    options_lat_lon_lines = {
        "stroke_width": 2.0,
        "layer": "meta",
        "stroke": COLOR_SCREWHOLES,
        "opacity": 0
    }

    options_lat_lon_lines_text = {
        "stroke_width": 0.5,
        "layer": "meta_text",
        "stroke": COLOR_SCREWHOLES,
        "opacity": 0
    }

    # options_text = {
    #     "stroke_width": 0.5,
    #     "layer": "meta_text",
    #     "stroke": [0, 0, 0],
    # }

    # if DARK_MODE:
    #     options["stroke"] = [255, 255, 255]
    #     options_text["stroke"] = [255, 255, 255]

    # options_tile = {
    #     "stroke_width": 4.0,
    #     "layer": "meta",
    #     "stroke": [0, 0, 0],
    #     "opacity": 0
    # }

    # -------------------------------------------------------------------------------------------------
    # draw screw holes

    if DRAW_SCREWHOLES:
        positions = [
            [MAP_FRAGMENT_OFFSET[0]+HOLE_DIST,                MAP_FRAGMENT_OFFSET[1]+HOLE_DIST],
            [MAP_FRAGMENT_OFFSET[0]+TILE_SIZE[0]-HOLE_DIST,   MAP_FRAGMENT_OFFSET[1]+HOLE_DIST],
            [MAP_FRAGMENT_OFFSET[0]+TILE_SIZE[0]-HOLE_DIST,   MAP_FRAGMENT_OFFSET[1]+TILE_SIZE[1]-HOLE_DIST],
            [MAP_FRAGMENT_OFFSET[0]+HOLE_DIST,                MAP_FRAGMENT_OFFSET[1]+TILE_SIZE[1]-HOLE_DIST],
        ]

        for col in range(0, TILE_NUMBERS[0]):
            for row in range(0, TILE_NUMBERS[1]):

                tile_origin = [TILE_SIZE[0]*col, TILE_SIZE[1]*row]

                for pos in positions:

                    x = tile_origin[0]+pos[0]
                    y = tile_origin[1]+pos[1]

                    p = Point([x, y])

                    if not viewport_polygon.contains(p):
                        continue

                    for i, _ in enumerate(tiles_flatten):
                        svg_handler.add_line([[x-2, y-2], [x+2, y+2]], **options_screwholes)
                        svg_handler.add_line([[x+2, y-2], [x-2, y+2]], **options_screwholes)

                        # svg.add_polygon(p.buffer(2), **options)
                        svg_handler.add_polygon(p.buffer(3+1), **options_screwholes)
                    exclusion_zones.append(p.buffer(7))

    # -------------------------------------------------------------------------------------------------
    # draw lat/lon lines
    # TODO: cut with viewport

    if DRAW_LAT_LON_LINES:
        latlonlines = []

        # color = [0, 0, 0]
        # if DARK_MODE:
        #     color = [255, 255, 255]

        NUM_LINES_LAT = 24*1

        for i in range(1, NUM_LINES_LAT): # lat 

            deg = 90 - (180/NUM_LINES_LAT)*i

            text_lines = get_text(hfont, "{:5.1f}".format(deg))
            text_lines = shapely.affinity.scale(text_lines, xfact=1, yfact=-1, origin=Point(0, 0))

            text_lines1 = shapely.affinity.translate(text_lines, xoff=2, yoff=MAP_SIZE[1]*i/NUM_LINES_LAT-1)
            for line in text_lines1.geoms:
                l = list(line.coords)

                if not viewport_polygon.contains(Point(l[0])) or not viewport_polygon.contains(Point(l[1])):
                    continue

                svg_handler.add_line(l, **options_lat_lon_lines_text)
            
            exclusion_zones.append(text_lines1.buffer(3).simplify(SIMPLIFICATION_MAX_ERROR))

            text_lines2 = shapely.affinity.translate(text_lines, xoff=MAP_FRAGMENT_SIZE[0]-20, yoff=MAP_SIZE[1]*i/NUM_LINES_LAT-1)
            for line in text_lines2.geoms:
                l = list(line.coords)

                if not viewport_polygon.contains(Point(l[0])) or not viewport_polygon.contains(Point(l[1])):
                    continue

                svg_handler.add_line(l, **options_lat_lon_lines_text)
            
            exclusion_zones.append(text_lines2.buffer(3).simplify(SIMPLIFICATION_MAX_ERROR))

            # exclusion_zones.append(Polygon())

            line = LineString([[0, MAP_SIZE[1]*i/NUM_LINES_LAT], [1+20, MAP_SIZE[1]*i/NUM_LINES_LAT]])
            exclusion_zones.append(line.buffer(2))
            latlonlines.append(line)

            line = LineString([[MAP_SIZE[0]-20, MAP_SIZE[1]*i/NUM_LINES_LAT], [MAP_SIZE[0], MAP_SIZE[1]*i/NUM_LINES_LAT]])
            exclusion_zones.append(line.buffer(2))
            latlonlines.append(line)

        NUM_LINES_LON = 24*2

        for i in range(1, NUM_LINES_LON): # lon

            line_length = 10

            if i % 2 == 0:
                line_length = 20

                deg = (180/NUM_LINES_LON)*i

                text_lines = get_text(hfont, "{:5.1f}".format(deg))
                text_lines = shapely.affinity.scale(text_lines, xfact=1, yfact=-1, origin=Point(0, 0))

                text_lines1 = shapely.affinity.translate(text_lines, xoff=MAP_SIZE[0]*i/NUM_LINES_LON+2, yoff=MAP_FRAGMENT_OFFSET[1]+20)
                for line in text_lines1.geoms:
                    l = list(line.coords)

                    if not viewport_polygon.contains(Point(l[0])) or not viewport_polygon.contains(Point(l[1])):
                        continue

                    svg_handler.add_line(l, **options_lat_lon_lines_text)

                exclusion_zones.append(text_lines1.buffer(3).simplify(SIMPLIFICATION_MAX_ERROR))

                text_lines2 = shapely.affinity.translate(text_lines, xoff=MAP_SIZE[0]*i/NUM_LINES_LON+2, yoff=MAP_FRAGMENT_OFFSET[1]+MAP_FRAGMENT_SIZE[1]-14)
                for line in text_lines2.geoms:
                    l = list(line.coords)

                    if not viewport_polygon.contains(Point(l[0])) or not viewport_polygon.contains(Point(l[1])):
                        continue

                    svg_handler.add_line(l, **options_lat_lon_lines_text)

                exclusion_zones.append(text_lines2.buffer(3).simplify(SIMPLIFICATION_MAX_ERROR))

            line = LineString([[MAP_SIZE[0]*i/NUM_LINES_LON, MAP_FRAGMENT_OFFSET[1]], [MAP_SIZE[0]*i/NUM_LINES_LON, MAP_FRAGMENT_OFFSET[1]+line_length]])
            exclusion_zones.append(line.buffer(2))
            latlonlines.append(line)

            line = LineString([[MAP_SIZE[0]*i/NUM_LINES_LON, MAP_FRAGMENT_OFFSET[1]+MAP_FRAGMENT_SIZE[1]-line_length], [MAP_SIZE[0]*i/NUM_LINES_LON, MAP_FRAGMENT_OFFSET[1]+MAP_FRAGMENT_SIZE[1]]])
            exclusion_zones.append(line.buffer(2))
            latlonlines.append(line)

        for line in latlonlines:

            save_line = line.intersection(viewport_polygon)
            if save_line.is_empty:
                continue

            svg_handler.add_line(save_line.coords, **options_lat_lon_lines)

    return

def draw_sea_labels(svg_handler):

    options_sea_labels = {
        "stroke_width": 1.5,
        "layer": "sea_labels",
        "stroke": COLOR_SEA_LABELS,
        "opacity": 0
    }

    for l in SEA_LABELS:

        pos, label = l

        if len(pos) == 0:
            continue

        c = conv.convert_wgs_to_map(*pos)

        text_lines = get_text(hfont_large, label)
        text_lines = shapely.affinity.scale(text_lines, xfact=1, yfact=-1, origin=Point(0, 0))
        text_lines = shapely.affinity.translate(text_lines, xoff=c[0], yoff=c[1])

        for line in text_lines.geoms:
            l = list(line.coords)
            if not viewport_polygon.contains(Point(l[0])) or not viewport_polygon.contains(Point(l[1])):
                continue

            color = [0, 0, 0]
            if DARK_MODE:
                color = [255, 255, 255]

            svg_handler.add_line(l, **options_sea_labels)
        
        exclusion_zones.append(text_lines.buffer(4+1).buffer(-1).simplify(SIMPLIFICATION_MAX_ERROR))

    return

def draw_coastlines(svg_handler):

    options_coastlines = {
        "stroke_width": 1.5,
        "layer": "coastlines",
        "stroke": COLOR_COASTLINES,
    }

    options_coastlines_hatching = {
        "stroke_width": 0,
        "layer": "coastlines_hatching",
        "hatching": "coastline_hatching"
    }

    coastlines = load_coastline()

    timer_start = datetime.now()

    coastlines_processed = []
    for i in range(0, len(coastlines)):
        c = coastlines[i]

        c = c.intersection(viewport_polygon)

        if c.area < 10:
            continue

        c = c.buffer(0.3).buffer(-0.3).simplify(SIMPLIFICATION_MAX_ERROR)

        if c.is_valid:
            coastlines_processed.append(c)
        else:
            print("error during postprocessing. coastline {}/{}".format(i, len(coastlines)))

    coastlines = coastlines_processed

    print(TIMER_STRING.format("postprocessing coastline data", (datetime.now()-timer_start).total_seconds()))

    for i, item in enumerate(tiles_flatten):
        tile = box(item[0][0], item[0][1], item[1][0], item[1][1])

        coastlines_combined = ops.unary_union(coastlines)
        coastlines_extended = coastlines_combined.buffer(4).difference(coastlines_combined)
        coastlines_extended = coastlines_extended.intersection(viewport_polygon)
        coastlines_extended = coastlines_extended.simplify(SIMPLIFICATION_MAX_ERROR)
        coastlines_extended = cut_polygons(coastlines_extended, exclusion_zones)

        for poly in coastlines_extended:
            if poly.area < 3:
                continue
            if type(poly) not in [Polygon, MultiPolygon]:
                print("warning: unknown geometry: {}".format(poly[0]))
                continue
            svg_handler.add_polygon(poly, **options_coastlines_hatching)

        coastlines_line = polygons_to_linestrings(coastlines)
        coastlines_line = cut_linestrings(coastlines_line, exclusion_zones)

        for coastline in coastlines_line:
            lines = []

            if coastline.is_empty:
                continue

            if type(coastline) is MultiLineString:
                for g in coastline.geoms:
                    lines.append(g)
            elif type(coastline) is LineString:
                lines.append(coastline)
            else:
                print("error: unknown geometry: {}".format(coastline[0]))

            for line in lines:

                color = [0, 0, 0]
                if DARK_MODE:
                    color = [255, 255, 255]

                svg_handler.add_poly_line(list(line.coords), **options_coastlines)

    return coastlines

def draw_bathymetry(svg_handler, cut_bathymetry_by):

    options = {
        "stroke_width": 0.2,
        "opacity": 0.0,
        "stroke": COLOR_BATHYMETRY,
        "layer": "bathymetry"
    }

    timer_start = datetime.now()

    bathymetry = []

    format_options = [BATHYMETRY_MIN_HEIGHT, BATHYMETRY_MAX_HEIGHT, BATHYMETRY_NUM_LAYERS]

    for i in range(0, BATHYMETRY_NUM_LAYERS):
        bathymetry.append([])

    for filename in BATHYMETRY_FILES:
        data = load_bathymetry_file(filename + "{}_{}_{}.geojson".format(*format_options), difference=cut_bathymetry_by)

        for i in range(0, len(data)):
            bathymetry[i] = bathymetry[i] + data[i]


    print(TIMER_STRING.format("parsing sea data", (datetime.now()-timer_start).total_seconds()))     
    timer_start = datetime.now()

    add_layers_to_writer(bathymetry, svg_handler, exclusion_zones)
    for item in [item for sublist in bathymetry for item in sublist]:
        exclusion_zones.append(item.simplify(SIMPLIFICATION_MAX_ERROR))

    print(TIMER_STRING.format("preparing sea data", (datetime.now()-timer_start).total_seconds()))

    return

def draw_terrain(svg_handler):

    options = {
        "stroke_width": 0.5,
        "stroke": COLOR_TERRAIN,
        "layer": "terrain"
    }

    timer_start = datetime.now()

    terrain = []

    format_options = [TERRAIN_MIN_HEIGHT, TERRAIN_MAX_HEIGHT, TERRAIN_NUM_LAYERS]

    for i in range(0, TERRAIN_NUM_LAYERS):
        terrain.append([])

    for filename in TERRAIN_FILES:
        data = load_bathymetry_file(filename + "{}_{}_{}.geojson".format(*format_options))

        for i in range(0, len(data)):
            terrain[i] = terrain[i] + data[i]

    print(TIMER_STRING.format("parsing terrain data", (datetime.now()-timer_start).total_seconds()))     
    timer_start = datetime.now()  

    # Tile border smoothing
    # for i in range(0, len(terrain)):

    #     # Expand every polygon by a bit so polygon merging (to remove tile-borders) works
    #     # at the same time this smoothes the rather rough lines a bit

    #     for j in range(0, len(terrain[i])):
    #         terrain[i][j] = terrain[i][j].buffer(+1.0).buffer(-0.8)

    #     terrain[i] = polygons_merge_tiles(terrain[i])

    # polygons to linestrings
    terrain = polygons_to_linestrings(terrain, flatten=False)

    print(TIMER_STRING.format("converting terrain data", (datetime.now()-timer_start).total_seconds()))  
    timer_start = datetime.now()  

    for i in range(0, len(terrain)):
        layer = terrain[i]
        for j in range(0, len(layer)):
            layer[j] = layer[j].intersection(viewport_polygon)



    print(TIMER_STRING.format("postprocessing terrain data", (datetime.now()-timer_start).total_seconds())) 

    for i in range(0, len(terrain)):
        layer = terrain[i]
        color = 0 # 0 + i * 3

        if DARK_MODE:
            color = 255

        layer_cut = cut_linestrings(layer, exclusion_zones)

        for terrain_line in layer_cut:
            lines = validate_linestring(terrain_line)

            for line in lines:
                svg_handler.add_poly_line(list(line.coords), **options)

    return

def draw_borders(svg_handler):

    options = {
        "stroke_width": 0.6,
        "stroke": COLOR_BORDERS,
        "layer": "borders"
    }

    timer_start = datetime.now()
    
    borders = get_lines_from_shapefile(BORDER_FILE, latlons_flipped=True)

    borders_cut = cut_linestrings(borders, exclusion_zones)

    for border_line in borders_cut:
        lines = validate_linestring(border_line)
        for line in lines:
            if viewport_polygon.contains(line):
                exclusion_zones.append(line.buffer(1).simplify(SIMPLIFICATION_MAX_ERROR))
                svg_handler.add_poly_line(list(line.coords), **options)

    print(TIMER_STRING.format("loading border region data", (datetime.now()-timer_start).total_seconds())) 

    return

def draw_cities(svg_handler):

    options = {
        "stroke_width": 0.5,
        "opacity": 0.0,
        "stroke": COLOR_CITIES_CIRCLES,
        "layer": "cities_circles"
    }

    options_text = {
        "stroke_width": 0.5,
        "opacity": 0.0,
        "stroke": COLOR_CITIES,
        "layer": "cities"
    }

    cities                   = []
    cities_names             = []
    city_name                = []
    city_pos                 = []
    city_label               = []
    cities_label_orientation = []

    geometries = []
    shapefile = fiona.open(CITIES_FILE)

    latlons_flipped = True

    if latlons_flipped:
        func = conv.convert_wgs_to_map_list_lon_lat
    else:
        func = conv.convert_wgs_to_map_list

    for item in shapefile:
        shapefile_geom = shape(item["geometry"])
        geom = ops.transform(func, shapefile_geom)
        geom = geom.simplify(SIMPLIFICATION_MAX_ERROR)

        if not type(geom) is Point:
            raise Exception("parsing shapefile: unexpected type: {}".format(geom))

        cities.append(geom)
        cities_names.append(item["properties"]["nameascii"])
        cities_label_orientation.append(item["properties"]["labelorientation"])


    for i in range(0, len(cities)):
        city_pos = cities[i]
        city_name = cities_names[i]
        city_label_orientation = cities_label_orientation[i]
        if viewport_polygon.contains(city_pos.buffer(CITY_CIRCLE_RADIUS)):
            svg_handler.add_polygon(city_pos.buffer(CITY_CIRCLE_RADIUS), **options)

        c = list(city_pos.coords)[0]

        text_lines = get_text(hfont, city_name)

        text_lines = shapely.affinity.scale(text_lines, xfact=1, yfact=-1, origin=Point(0, 0))

        minx, miny, maxx, maxy = text_lines.bounds
        text_width = maxx-minx
        text_height = maxy-miny

        if city_label_orientation is None:
            city_label_orientation = "ne"
        if city_label_orientation == "n":
            x_offset = c[0]-(text_width/2)
            y_offset = c[1]-2
        elif city_label_orientation == "ne":
            x_offset = c[0]+CITY_CIRCLE_RADIUS*2
            y_offset = c[1]-1
        elif city_label_orientation == "e":
            x_offset = c[0]+CITY_CIRCLE_RADIUS*2.5
            y_offset = c[1]
        elif city_label_orientation == "se":
            x_offset = c[0]+CITY_CIRCLE_RADIUS*2
            y_offset = c[1]+1+text_height
        elif city_label_orientation == "s":
            x_offset = c[0]-(text_width/2)
            y_offset = c[1]+2+text_height
        elif city_label_orientation == "sw":
            x_offset = c[0]-(text_width)-(CITY_CIRCLE_RADIUS*2)
            y_offset = c[1]+1+text_height
        elif city_label_orientation == "w":
            x_offset = c[0]-(text_width)-(CITY_CIRCLE_RADIUS*2.5)
            y_offset = c[1]+(text_height/2)
        elif city_label_orientation == "nw":
            x_offset = c[0]-(text_width)-(CITY_CIRCLE_RADIUS*2)
            y_offset = c[1]-1

        # text_lines = shapely.affinity.translate(text_lines, xoff=c[0]+CITY_CIRCLE_RADIUS*2, yoff=c[1]-1.0)
        text_lines = shapely.affinity.translate(text_lines, xoff=x_offset, yoff=y_offset)

        for line in text_lines.geoms:
            if viewport_polygon.contains(line):
                l = list(line.coords)
                svg_handler.add_line(l, **options_text)

        exclusion_zones.append(text_lines.buffer(2+1).buffer(-1).simplify(SIMPLIFICATION_MAX_ERROR))
        exclusion_zones.append(city_pos.buffer(CITY_CIRCLE_RADIUS+1).simplify(SIMPLIFICATION_MAX_ERROR))

    return

def create_cache():
    timer_start = datetime.now()

    print("writing bathymetric data to cache")

    coastlines = []

    viewport_polygon = Polygon([
        [VIEWPORT_OFFSET[0],                    VIEWPORT_OFFSET[1]],
        [VIEWPORT_OFFSET[0]+VIEWPORT_SIZE[0],   VIEWPORT_OFFSET[1]],
        [VIEWPORT_OFFSET[0]+VIEWPORT_SIZE[0],   VIEWPORT_OFFSET[1]+VIEWPORT_SIZE[1]],
        [VIEWPORT_OFFSET[0],                    VIEWPORT_OFFSET[1]+VIEWPORT_SIZE[1]]
    ])

    coastlines = load_coastline()
    timer_start = datetime.now()
    coastlines_processed = []
    for i in range(0, len(coastlines)):
        c = coastlines[i]
        c = c.intersection(viewport_polygon)
        if c.area < 10:
            continue
        c = c.buffer(0.3).buffer(-0.3).simplify(SIMPLIFICATION_MAX_ERROR)
        if c.is_valid:
            coastlines_processed.append(c)
        else:
            print("error during postprocessing. coastline {}/{}".format(i, len(coastlines)))
    coastlines = coastlines_processed
    # num_layers = 15
    # min_height = -9000
    # max_height = 0
    format_options = [BATHYMETRY_MIN_HEIGHT, BATHYMETRY_MAX_HEIGHT, BATHYMETRY_NUM_LAYERS]
    for filename in BATHYMETRY_FILES:
        write_bathymetry_to_cache(filename + "{}_{}_{}.geojson".format(*format_options), coastlines)


    # TODO: remove necessity to indicate num_layer, min_height and max_height
    # terrain
    # num_layers = 30
    # min_height = 0
    # max_height = 9000

    format_options = [TERRAIN_MIN_HEIGHT, TERRAIN_MAX_HEIGHT, TERRAIN_NUM_LAYERS]
    for filename in BATHYMETRY_FILES:
        write_bathymetry_to_cache(filename + "{}_{}_{}.geojson".format(*format_options))

    print(TIMER_STRING.format("creating cache", (datetime.now()-timer_start).total_seconds())) 

    return

# --------------------------------------------------------------------------------
# MAIN ---->

timer_total = datetime.now()

conv = maptools.Converter(MAP_CENTER, MAP_SIZE, MAP_SIZE_SCALE)

# --------------------------------------------------------------------------------
# create flattened tile matrix e.g. 3x2 --> [[tile-1-1, tile-1-2, tile-1-3, tile-2-1, tile-2-2, tile-2-3]]
# tile-row-column = [(minx, miny), (maxx, maxy)]
x_range = [TILE_SIZE[0]*i for i in range(TILE_NUMBERS[0]+1)]
y_range = [TILE_SIZE[1]*i for i in range(TILE_NUMBERS[1]+1)]
tiles =[]
for i, item_y in enumerate(y_range[:-1]):
    row = []
    for j, item_x in enumerate(x_range[:-1]):
        col = [(x_range[j], y_range[i]), (x_range[j+1], y_range[i+1])]
        row.append(col)
    tiles.append(row)
tiles_flatten = [item for sublist in tiles for item in sublist]
print(tiles_flatten)

# --------------------------------------------------------------------------------
# fonts
hfont = HersheyFonts()
hfont.load_default_font("futural")
hfont.normalize_rendering(FONT_SIZE)
hfont_large = HersheyFonts()
hfont_large.load_default_font("futural")
hfont_large.normalize_rendering(FONT_SIZE_LARGE)

# --------------------------------------------------------------------------------
# create cache of the complete viewport
# TODO: create cache automatically when no files found
if CREATE_CACHE:
    create_cache()

# --------------------------------------------------------------------------------
# instead of creating svg(s) of the whole map, create a small svg-preview
if PREVIEW:

    print("map size: {:.2f} x {:.2f} meter".format(MAP_SIZE[0]*MAP_SIZE_SCALE, MAP_SIZE[1]*MAP_SIZE_SCALE))
    print("svg size: {:.2f} x {:.2f} units".format(*MAP_SIZE))
    print("viewport size: {} x {} millimeter".format(*VIEWPORT_SIZE))
    print("creating a preview")
    print("preview size: {} x {} millimeter".format(*PREVIEW_SIZE))
    print("preview position: {} x {} millimeter".format(*PREVIEW_POSITION))

    svg_preview = SvgWriter("world-preview.svg", dimensions=PREVIEW_SIZE, offset=[PREVIEW_POSITION[0]+VIEWPORT_OFFSET[0], PREVIEW_POSITION[1]+VIEWPORT_OFFSET[1]], background_color=BG_COLOR)

    exclusion_zones = []
    coastlines = []
    places = []
    urban = []
    borders = []
    admin = []
    bathymetry = []
    terrain = []

    svg_preview.add_layer("urban")
    svg_preview.add_layer("borders")
    svg_preview.add_layer("bathymetry")
    svg_preview.add_layer("terrain")
    svg_preview.add_layer("coastlines")
    svg_preview.add_layer("coastlines_hatching")
    svg_preview.add_layer("places")
    svg_preview.add_layer("places_circles")
    svg_preview.add_layer("cities")
    svg_preview.add_layer("cities_circles")
    svg_preview.add_layer("sea_labels")
    svg_preview.add_layer("meta")
    svg_preview.add_layer("meta_text")
    svg_preview.add_hatching("coastline_hatching", stroke_width=0.5, distance=2.0)

    # -------------------------------------------------------------------------------------------------------
    # create bathymetry hatching
    for i in range(0, BATHYMETRY_NUM_LAYERS):

        if i < 5:
            distance = 1
        else:
            distance = 1.15**(i-4)

        svg_preview.add_hatching("bathymetry_hatching_{}".format(i), stroke_width=0.5, stroke_opacity=0.5, distance=distance)
    
    # -------------------------------------------------------------------------------------------------------
    viewport_polygon = Polygon([
        [VIEWPORT_OFFSET[0]+PREVIEW_POSITION[0],                VIEWPORT_OFFSET[1]+PREVIEW_POSITION[1]],
        [VIEWPORT_OFFSET[0]+PREVIEW_POSITION[0]+PREVIEW_SIZE[0], VIEWPORT_OFFSET[1]+PREVIEW_POSITION[1]],
        [VIEWPORT_OFFSET[0]+PREVIEW_POSITION[0]+PREVIEW_SIZE[0], VIEWPORT_OFFSET[1]+PREVIEW_POSITION[1]+PREVIEW_SIZE[1]],
        [VIEWPORT_OFFSET[0]+PREVIEW_POSITION[0],                VIEWPORT_OFFSET[1]+PREVIEW_POSITION[1]+PREVIEW_SIZE[1]]
    ])

    # -------------------------------------------------------------------------------------------------------
    if DRAW_META:
        draw_meta(svg_preview)
    if DRAW_SEA_LABELS:
        draw_sea_labels(svg_preview)
    if DRAW_CITIES:
        draw_cities(svg_preview)
    if DRAW_COASTLINES:
        coastlines = draw_coastlines(svg_preview)
    if DRAW_BORDERS:
        draw_borders(svg_preview)
    if DRAW_BATHYMETRY:
        draw_bathymetry(svg_preview, coastlines)
    if DRAW_TERRAIN:
        draw_terrain(svg_preview)

    svg_preview.save() 
    print("total time   ->   {}s".format((datetime.now()-timer_total).total_seconds()))
    sys.exit()

# --------------------------------------------------------------------------------
# create svg_handler for each item in tiles
svg = []
for i, item in enumerate(tiles_flatten):
    svg.append(SvgWriter("world-{}.svg".format(i), dimensions=TILE_SIZE, offset=[item[0][0]+VIEWPORT_OFFSET[0], item[0][1]+VIEWPORT_OFFSET[1]], background_color=BG_COLOR))

# --------------------------------------------------------------------------------
print("map size: {:.2f} x {:.2f} meter".format(MAP_SIZE[0]*MAP_SIZE_SCALE, MAP_SIZE[1]*MAP_SIZE_SCALE))
print("svg size: {:.2f} x {:.2f} units".format(*MAP_SIZE))
print("viewport size: {:.2f} x {:.2f} millimeter".format(*VIEWPORT_SIZE))
print("tile numbers: {} x {}".format(*TILE_NUMBERS))
print("tile size: {:.2f} x {:.2f} millimeter".format(*TILE_SIZE))

# --------------------------------------------------------------------------------
# loop through each tile, create the map and save as svg
for tile_index, tile_item in enumerate(tiles_flatten):

    timer_tile = datetime.now()
    print("starting tile {} of {}".format(tile_index + 1, len(tiles_flatten)))

    exclusion_zones = []
    coastlines = []
    places = []
    urban = []
    borders = []
    admin = []
    bathymetry = []
    terrain = []

    svg[tile_index].add_layer("urban")
    svg[tile_index].add_layer("borders")
    svg[tile_index].add_layer("bathymetry")
    svg[tile_index].add_layer("terrain")
    svg[tile_index].add_layer("coastlines")
    svg[tile_index].add_layer("coastlines_hatching")
    svg[tile_index].add_layer("places")
    svg[tile_index].add_layer("places_circles")
    svg[tile_index].add_layer("cities")
    svg[tile_index].add_layer("cities_circles")
    svg[tile_index].add_layer("sea_labels")
    svg[tile_index].add_layer("meta")
    svg[tile_index].add_layer("meta_text")
    svg[tile_index].add_hatching("coastline_hatching", stroke_width=0.5, distance=2.0)

    # -------------------------------------------------------------------------------------------------------
    # create bathymetry hatching
    if DARK_MODE:
        for i in range(0, BATHYMETRY_NUM_LAYERS):

            if i < 5:
                distance = 6
            else:
                distance = 1.5 + 1.2**(14-4) - (1.2**(i-4) - 1.2**(5-4))

            svg[tile_index].add_hatching("bathymetry_hatching_{}".format(i), stroke_width=0.5, stroke_opacity=0.5, distance=distance)
    else:
        for i in range(0, BATHYMETRY_NUM_LAYERS):

            if i < 5:
                distance = 1
            else:
                distance = 1.15**(i-4)

            svg[tile_index].add_hatching("bathymetry_hatching_{}".format(i), stroke_width=0.5, stroke_opacity=0.5, distance=distance)
    # -------------------------------------------------------------------------------------------------------

    viewport_polygon = Polygon([
        [VIEWPORT_OFFSET[0]+tile_item[0][0],   VIEWPORT_OFFSET[1]+tile_item[0][1]],
        [VIEWPORT_OFFSET[0]+tile_item[1][0],   VIEWPORT_OFFSET[1]+tile_item[0][1]],
        [VIEWPORT_OFFSET[0]+tile_item[1][0],   VIEWPORT_OFFSET[1]+tile_item[1][1]],
        [VIEWPORT_OFFSET[0]+tile_item[0][0],   VIEWPORT_OFFSET[1]+tile_item[1][1]]
    ])

    if DRAW_META:
        draw_meta(svg[tile_index])
    if DRAW_SEA_LABELS:
        draw_sea_labels(svg[tile_index])
    if DRAW_CITIES:
        draw_cities(svg[tile_index])
    if DRAW_COASTLINES:
        coastlines = draw_coastlines(svg[tile_index])
    if DRAW_BORDERS:
        draw_borders(svg[tile_index])
    if DRAW_BATHYMETRY:
        draw_bathymetry(svg[tile_index], coastlines)
    if DRAW_TERRAIN:
        draw_terrain(svg[tile_index])

    svg[tile_index].save() 
    print("time for tile {} of {}   ->   {}s".format(tile_index + 1, len(tiles_flatten), (datetime.now()-timer_tile).total_seconds()))

print("total time   ->   {}s".format((datetime.now()-timer_total).total_seconds()))

# <---- MAIN
# --------------------------------------------------------------------------------

# if DRAW_URBAN_AREAS:

#     URBAN_FILE = "world_data/10m_cultural/ne_10m_urban_areas.shp"

#     timer_start = datetime.now()
    
#     urban = get_polys_from_shapefile(URBAN_FILE, latlons_flipped=True)

#     print(TIMER_STRING.format("loading urban areas data", (datetime.now()-timer_start).total_seconds())) 

#     timer_start = datetime.now()

#     areas_postprocessed = []
#     for i in range(0, len(urban)):
#         a = urban[i]
#         a = a.buffer(0.1)

#         # if a.area < 2:
#         #     continue

#         areas_postprocessed.append(a)

#     urban = []
#     unified_poly = ops.unary_union(areas_postprocessed)
#     for poly in unified_poly.geoms:
#         if poly.area > 8:
#             urban.append(poly.buffer(1.0).buffer(1.0))

#     print(TIMER_STRING.format("postprocessing urban areas data", (datetime.now()-timer_start).total_seconds())) 

# --------------------------------------------------------------------------------

# if DRAW_ADMIN_REGIONS:

#     ADMIN_REGIONS_FILE = "world_data/10m_cultural/ne_10m_admin_1_states_provinces_lines.shp"

#     timer_start = datetime.now()
    
#     admin = get_lines_from_shapefile(ADMIN_REGIONS_FILE, latlons_flipped=True)

#     print(TIMER_STRING.format("loading admin region data", (datetime.now()-timer_start).total_seconds())) 


# --------------------------------------------------------------------------------

# for place in places:
#     svg.add_circles(list(place.coords), radius=0.5, layer="places")

# ---

# for u in urban:
#     if not u.is_valid:
#         continue
#     svg.add_polygon(u, stroke_width=0.2, opacity=0.3, stroke=[255, 255, 255], layer="urban")

# ---

# for a in admin:
#     svg.add_poly_line(list(a.coords), stroke_width=0.2, stroke=[255, 255, 255], layer="borders")

