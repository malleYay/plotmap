import math
from datetime import datetime
import random

from shapely.geometry import GeometryCollection, LineString, MultiLineString, Polygon, MultiPolygon, Point, MultiPoint
from shapely.ops import unary_union

class SvgWriter(object):

    def __init__(self, filename, dimensions=None, offset=None, image=None, background_color=None):

        if not filename.endswith(".svg"):
            filename += ".svg"

        self.filename   = filename
        self.dimensions = dimensions
        self.offset     = offset
        if self.offset is None:
            self.offset = [0, 0]
        self.image      = image
        self.background_color = background_color

        self.hatchings                  = {} # hatchings MultiLineStrings
        self.hatching_options           = {} # hatchings SVG options (stroke-width, ...)
        self.hatching_options_meta      = {} # hatchings meta options (distance, orientation, ...)

        self.layers     = {}
        self.add_layer("default")

    def add_layer(self, layer_id):
        self.layers[layer_id] = {}
        self.layers[layer_id]["hexagons"]   = []
        self.layers[layer_id]["circles"]    = []
        self.layers[layer_id]["rectangles"] = []
        self.layers[layer_id]["polygons"]   = []
        self.layers[layer_id]["lines"]      = []
        self.layers[layer_id]["poly_lines"] = []
        self.layers[layer_id]["raw"]        = []

    def add_circles(self, circles, radius=3, fill=[255, 0, 0], layer="default"):
        for item in circles:
            self.layers[layer]["circles"].append([item, radius, fill])

    # coords: [[x, y], [width, height]]
    def add_rectangle(self, coords, stroke_width=1, stroke=[255, 0, 0], opacity=1.0, layer="default"):
        self.layers[layer]["rectangles"].append([*coords, stroke_width, stroke, opacity])

    # coords: [[x1, y1], [x2, x2]]
    def add_line(self, coords, stroke_width=1, stroke=[0, 0, 0], stroke_opacity=1.0, stroke_dasharray=None, opacity=None, layer="default"):

        if len(coords) != 2:
            # raise Exception("add_line: malformed input data: {}".format(coords))
            raise Exception("add_line: malformed input data")

        options = {}
        options["stroke-width"]     = stroke_width
        options["stroke"]           = stroke
        options["stroke-opacity"]   = stroke_opacity

        # parameter opacity is ignored and only exists to ensure compatibility with the polygon functions
        # (so the same options hashmap can be used on all of these functions)

        if stroke_dasharray is not None:
            options["stroke-dasharray"] = stroke_dasharray

        self.layers[layer]["lines"].append((coords, options))

    def add_lines(self, coords, **kwargs):
        for coord in coords:
            self.add_line(coord, **kwargs)

    # poly may be Polygon or MultiPolygon
    def add_polygon(self, poly, 
            stroke_width=1, 
            stroke=[0, 0, 0], 
            fill=[120, 120, 120], 
            opacity=1.0, 
            repeat=1,
            layer="default", 
            hatching=None):

        options = {}
        options["stroke-width"]     = stroke_width
        options["stroke"]           = stroke
        options["fill"]             = fill
        options["opacity"]          = opacity

        polys = []

        if type(poly) is Polygon:
            polys.append(poly)
        elif type(poly) is MultiPolygon:
            polys += list(poly.geoms)
        elif type(poly) is GeometryCollection:
            for g in poly.geoms:
                if type(g) is Polygon:
                    polys.append(poly)
                elif type(g) is MultiPolygon:
                    polys += list(poly.geoms)
        else:
            # raise Exception("unknown geometry: {}".format(poly))
            raise Exception("unknown geometry - add polygon - type {}".format(type(poly)))

        for p in polys:
            if stroke_width > 0:
                self.layers[layer]["polygons"].append((p.exterior.coords, options))
                for hole in p.interiors:
                    self.layers[layer]["polygons"].append((hole.coords, options))

            if hatching is not None:
                kwargs = {}
                kwargs["stroke_width"]  = stroke_width
                kwargs["stroke"]        = stroke
                kwargs["layer"]         = layer
                self._add_hatching_for_polygon(p, hatching, kwargs)

    def add_poly_line(self, coords, stroke_width=1, stroke=[0, 0, 0], stroke_opacity=1.0, opacity=None, layer="default"):
        options = {}
        options["stroke-width"]     = stroke_width
        options["stroke"]           = stroke
        options["stroke-opacity"]   = stroke_opacity

        # parameter opacity is ignored and only exists to ensure compatibility with the polygon functions
        # (so the same options hashmap can be used on all of these functions)

        self.layers[layer]["poly_lines"].append((coords, options))

    def add_raw_element(self, text, layer="default"):
        self.layers[layer]["raw"].append(text)

    @staticmethod
    def _line_intersection(line1, line2): 

        A = line1[0]
        B = line1[1]
        C = line2[0]
        D = line2[1]

        Bx_Ax = B[0] - A[0] 
        By_Ay = B[1] - A[1] 
        Dx_Cx = D[0] - C[0] 
        Dy_Cy = D[1] - C[1] 
        
        determinant = (-Dx_Cx * By_Ay + Bx_Ax * Dy_Cy) 
        
        if abs(determinant) < 1e-20: 
            return None 

        s = (-By_Ay * (A[0] - C[0]) + Bx_Ax * (A[1] - C[1])) / determinant 
        t = ( Dx_Cx * (A[1] - C[1]) - Dy_Cy * (A[0] - C[0])) / determinant 

        if s >= 0 and s <= 1 and t >= 0 and t <= 1: 
            return (A[0] + (t * Bx_Ax), A[1] + (t * By_Ay)) 

        return None

    # rotation: 45-90 degrees
    # TODO: 0-45 degrees

    HATCHING_ORIENTATION_45         = 0x01
    HATCHING_ORIENTATION_45_REV     = 0x02
    HATCHING_ORIENTATION_VERTICAL   = 0x03
    HATCHING_ORIENTATION_HORIZONTAL = 0x04

    def add_hatching(self, name, orientation=HATCHING_ORIENTATION_45, distance=2, bounding_box=None, stroke_width=0.2, stroke_dasharray=None, stroke_opacity=1.0, wiggle=0):

        self.hatching_options[name] = {}
        self.hatching_options[name]["stroke_width"] = stroke_width
        self.hatching_options[name]["stroke_dasharray"] = stroke_dasharray
        self.hatching_options[name]["stroke_opacity"] = stroke_opacity

        self.hatching_options_meta[name] = {}
        self.hatching_options_meta[name]["distance"] = distance
        self.hatching_options_meta[name]["orientation"] = orientation
        self.hatching_options_meta[name]["wiggle"] = wiggle

        # new_hatchlines = self._add_hatching(orientation=orientation, distance=distance, wiggle=wiggle, bounding_box=bounding_box)

        # if not name in self.hatchings:
        #     self.hatchings[name] = new_hatchlines
        # else:
        #     all_hatchlines = []

        #     for g in new_hatchlines.geoms:
        #         all_hatchlines.append(g)
        #     for g in self.hatchings[name]:
        #         all_hatchlines.append(g)

        #     self.hatchings[name] = MultiLineString(all_hatchlines)

    def _add_hatching(self, orientation=HATCHING_ORIENTATION_45, distance=2, wiggle=0, bounding_box=None):

        minx, miny, maxx, maxy = 0, 0, self.dimensions[0], self.dimensions[1]

        if bounding_box is not None:
            minx, miny, maxx, maxy = bounding_box
            minx = int(minx/distance) * distance
            miny = int(miny/distance) * distance
            maxx = (int(maxx/distance) + 1) * distance
            maxy = (int(maxy/distance) + 1) * distance

        height = maxy-miny
        width = maxx-minx

        num_lines = (width + height)/float(distance)

        if orientation == self.HATCHING_ORIENTATION_HORIZONTAL:
            num_lines = width/float(distance)

        if orientation == self.HATCHING_ORIENTATION_VERTICAL:
            num_lines = height/float(distance)

        north = [[minx, miny], [maxx, miny]]
        south = [[minx, maxy], [maxx, maxy]]
        west  = [[minx, miny], [minx, maxy]]
        east  = [[maxx, miny], [maxx, maxy]]

        hatchlines = []

        wiggle_range = [-wiggle, +wiggle]

        for i in range(0, int(num_lines)):

            random_error_1 = 0
            random_error_2 = 0

            if wiggle > 0:
                random_error_1 = random.uniform(*wiggle_range)
                random_error_2 = random.uniform(*wiggle_range)

            if orientation == self.HATCHING_ORIENTATION_45:
                x1 = minx
                y1 = miny + i * distance
                x2 = minx + i * distance
                y2 = miny
                y1 += random_error_1
                x2 += random_error_2
            elif orientation == self.HATCHING_ORIENTATION_45_REV:
                x1 = maxx
                y1 = miny + i * distance
                x2 = maxx - i * distance
                y2 = miny
                y1 += random_error_1
                x2 += random_error_2
            elif orientation == self.HATCHING_ORIENTATION_VERTICAL:
                x1 = minx + i * distance
                y1 = miny
                x2 = x1
                y2 = maxy
                x1 += random_error_1
                x2 += random_error_2
            elif orientation == self.HATCHING_ORIENTATION_HORIZONTAL:
                x1 = minx
                y1 = miny + i * distance
                x2 = maxx
                y2 = y1
                y1 += random_error_1
                y2 += random_error_2

            else:
                raise Exception("unknown hatching orientation type: {}".format(orientation))

            hatching_line = [[x1, y1], [x2, y2]]
            cropped_line = []

            north_intersect = SvgWriter._line_intersection(hatching_line, north)
            south_intersect = SvgWriter._line_intersection(hatching_line, south)
            west_intersect = SvgWriter._line_intersection(hatching_line, west)
            east_intersect = SvgWriter._line_intersection(hatching_line, east)

            if west_intersect is not None:
                cropped_line.append(west_intersect)

            if south_intersect is not None:
                cropped_line.append(south_intersect)

            if north_intersect is not None:
                cropped_line.append(north_intersect)

            if east_intersect is not None:
                cropped_line.append(east_intersect)

            if len(cropped_line) == 2:
                hatchlines.append(LineString(cropped_line))
            elif len(cropped_line) > 2:
                hatchlines.append(LineString([cropped_line[0], cropped_line[2]]))

        # if len(hatchlines) == 0:
        #     raise Exception("no hatchlines created for distance {}".format(distance))

        return MultiLineString(hatchlines)

    def _add_hatching_for_polygon(self, poly, hatching_name, polygon_options):

        hatchlines = []

        if (len(self.hatching_options_meta[hatching_name])) <= 0:
            # raise Exception("missing hatching: {}".format(hatching_name))
            raise Exception("missing hatching")
        
        # hatchlines = self.hatchings[hatching_name]
        hatchlines = self._add_hatching(**self.hatching_options_meta[hatching_name], bounding_box=poly.bounds)

        intersections = poly.intersection(hatchlines)

        if intersections.is_empty:
            return

        options = {
            **polygon_options,
            **self.hatching_options[hatching_name]
        }

        if type(intersections) is LineString:
            self.add_line(intersections.coords, **options)
        elif type(intersections) is MultiLineString:
            for line in intersections.geoms:
                self.add_line(line.coords, **options)
        elif type(intersections) is GeometryCollection:
            for line in intersections.geoms:
                if type(line) is LineString:
                    self.add_line(line.coords, **options)
                else:
                    pass
                    # print("unknown sub-geometry: {}".format(line))
        elif type(intersections) is MultiPoint:
            return
        elif type(intersections) is Point:
            return
        else:
            # raise Exception("error: unknown geometry: {}".format(type(intersections)))
            raise Exception("error: unknown geometry - add hatching")


    def save(self):

        timer_start = datetime.now()

        with open(self.filename, "w") as out:

            out.write("<?xml version=\"1.0\" encoding=\"utf-8\" ?>")
            out.write("<?xml-stylesheet href=\"style.css\" type=\"text/css\" title=\"main_stylesheet\" alternate=\"no\" media=\"screen\" ?>")
            
            if self.dimensions is not None:
                out.write("<svg baseProfile=\"tiny\" version=\"1.2\" width=\"{}px\" height=\"{}px\" ".format(self.dimensions[0], self.dimensions[1]))
            else:
                out.write("<svg baseProfile=\"tiny\" version=\"1.2\" ")
            out.write("xmlns=\"http://www.w3.org/2000/svg\" ")
            out.write("xmlns:ev=\"http://www.w3.org/2001/xml-events\" ")
            out.write("xmlns:xlink=\"http://www.w3.org/1999/xlink\" ")
            out.write("xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\" ")
            out.write(">")
            out.write("<defs />")

            if self.image is not None:
                out.write("<image x=\"0\" y=\"0\" xlink:href=\"{}\" />".format(self.image))

            if self.background_color is not None:
                out.write("<style>svg {{ background-color: {}; }}</style>".format(self.background_color))

            for layerid in self.layers.keys():

                # if layerid == "default":
                #     continue

                layer = self.layers[layerid]
                
                out.write("<g inkscape:groupmode=\"layer\" id=\"{0}\" inkscape:label=\"{0}\">".format(layerid))

                for c in layer["circles"]:
                    out.write("<circle cx=\"{}\" cy=\"{}\" fill=\"rgb({},{},{})\" r=\"{}\" />".format(
                        c[0][0]-self.offset[0], 
                        c[0][1]-self.offset[1], 
                        c[2][0], c[2][1], c[2][2], 
                        c[1]))

                for r in layer["rectangles"]:
                    out.write("<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" stroke-width=\"{}\" stroke=\"rgb({},{},{})\" fill-opacity=\"0.0\" stroke-opacity=\"{}\" />".format(
                        r[0][1]--self.offset[0], r[0][1]-self.offset[1],
                        *r[1], r[2], *r[3], r[4]))

                for line in layer["lines"]:
                    l = line[0]
                    options = line[1]
                    out.write("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" ".format(
                        l[0][0]-self.offset[0], l[0][1]-self.offset[1],
                        l[1][0]-self.offset[0], l[1][1]-self.offset[1]))
                    out.write("stroke-width=\"{}\" ".format(options["stroke-width"]))
                    out.write("stroke=\"rgb({}, {}, {})\" ".format(*options["stroke"]))

                    if "stroke-dasharray" in options:
                        out.write("stroke-dasharray=\"{}\" ".format(options["stroke-dasharray"]))

                    out.write("/>")

                for poly in layer["polygons"]:
                    p = poly[0]
                    options = poly[1]
                    out.write("<path d=\"")
                    out.write("M{} {} ".format(float(p[0][0]-self.offset[0]), float(p[0][1]-self.offset[1])))
                    for point in p[1:]:
                        out.write("L")
                        out.write(str(float(point[0]-self.offset[0])))
                        out.write(" ")
                        out.write(str(float(point[1]-self.offset[1])))
                        out.write(" ")
                    out.write("Z\" ")
                    out.write("stroke-width=\"{}\" ".format(options["stroke-width"]))
                    out.write("stroke=\"rgb({},{},{})\" ".format(*options["stroke"]))
                    out.write("fill=\"rgb({},{},{})\" ".format(*options["fill"]))
                    out.write("fill-opacity=\"{}\" />".format(options["opacity"]))

                for line in layer["poly_lines"]:
                    l = line[0]
                    options = line[1]
                    out.write("<path d=\"")
                    out.write("M{} {} ".format(float(l[0][0]-self.offset[0]), float(l[0][1]-self.offset[1])))
                    for point in l[1:]:
                        out.write("L")
                        out.write(str(float(point[0]-self.offset[0])))
                        out.write(" ")
                        out.write(str(float(point[1]-self.offset[1])))
                        out.write(" ")
                    out.write("\" ")
                    out.write("stroke-width=\"{}\" ".format(options["stroke-width"]))
                    out.write("stroke=\"rgb({},{},{})\" ".format(*options["stroke"]))
                    out.write("stroke-opacity=\"{}\" ".format(options["stroke-opacity"]))
                    out.write("fill=\"rgb({},{},{})\" ".format(0, 0, 0))
                    out.write("fill-opacity=\"{}\" />".format(0))
                    out.write("/>")

                for r in layer["raw"]:
                    out.write(r)

                out.write("</g>")

            out.write("</svg>")

        print("writing SVG in {0:.2f}s".format((datetime.now()-timer_start).total_seconds()))