from h3 import h3
import geopandas as gpd
from shapely.geometry import Polygon, Point
import numpy as np

def create_hexagon_around_point(h3_resolution, point):
    """
    Creates a hexagon around a point.

    Args:
        h3_resolution (int): H3 resolution level
        point (shapely.geometry.Point): Point around which the hexagon is to be created

    Returns:
        shapely.geometry.Polygon: Hexagon around the point
    """
    hexagon_vertices = h3.h3_to_geo_boundary(h3.geo_to_h3(point.y, point.x, h3_resolution))
    hexagon_polygon = Polygon(hexagon_vertices)
    return hexagon_polygon

def assign_hex_ids(data, h3_resolution):
    """
    Assigns H3 hexagon IDs to data points and calculates distances from centroids.

    Args:
        data (gpd.GeoDataFrame): GeoDataFrame containing the data points
        h3_resolution (int): H3 resolution level

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with hexagon IDs and distances from centroids
    """
    hexagons = []
    for _, row in data.iterrows():
        hexagon_polygon = create_hexagon_around_point(h3_resolution, row['geometry'],)
        hex_ids = h3.polyfill({'type': 'Polygon', 'coordinates': [list(hexagon_polygon.exterior.coords)]}, h3_resolution)
        hexagons.extend(hex_ids)

    data['hex_id'] = hexagons
    data['hex_centroid'] = data['hex_id'].apply(lambda x: Point(h3.h3_to_geo(x)))
    # data['distance_from_centroid'] = data.apply(lambda x: x['geometry'].distance(x['hex_centroid']), axis=1)
    hdf = gpd.GeoDataFrame(data['hex_id'], geometry=[Polygon(h3.h3_to_geo_boundary(h, geo_json=True)) for h in hexagons])
    return data, hdf

def containment_filter(h3_cell, points):
    """ Filter points within a hexagon

    Args:
        h3_cell (str): H3 cell ID
        points (gpd.GeoDataFrame): GeoDataFrame containing the points to filter

    Returns:
        gpd.GeoDataFrame: GeoDataFrame containing the points within the hexagon
    """
    boundary_vertices = h3.h3_to_geo_boundary(h3_cell, True)
    vlist = list(map(Point, boundary_vertices))
    bdf = gpd.GeoDataFrame(geometry=vlist, crs="EPSG:4326")
    polygon = Polygon(boundary_vertices)
    points = points.to_crs("EPSG:4326")
    return points[points.geometry.apply(lambda x: x.within(polygon))]

def containment_filter_generic(bdry, points):
    """ Filter points within a polygon

    Args:
        bdry (gpd.GeoDataFrame): GeoDataFrame containing the boundary polygon
        points (gpd.GeoDataFrame): GeoDataFrame containing the points to filter

    Returns:
        gpd.GeoDataFrame: GeoDataFrame containing the points within the polygon
    """
    polygon = bdry.geometry[0]
    points = points.to_crs(bdry.crs)
    return points[points.geometry.apply(lambda x: x.within(polygon))]

def generate_grid_points(h3_cell, epsg, grid_spacing):
    """ Generate grid points within a hexagon

    Args:
        h3_cell (str): H3 cell ID
        epsg (int): EPSG code for the coordinate reference system
        grid_spacing (float): Spacing between grid points

    Returns:
        np.array: x-coordinates of grid points
        np.array: y-coordinates of grid points
    """
    boundary_vertices = h3.h3_to_geo_boundary(h3_cell, True)
    vlist = list(map(Point, boundary_vertices))
    bdf = gpd.GeoDataFrame(geometry=vlist, crs="EPSG:4326")
    bdf = bdf.to_crs(epsg)
    converted_vertices = list(bdf.geometry)

    # Get bounding box
    min_lon = min([v.x for v in converted_vertices])
    max_lon = max([v.x for v in converted_vertices])
    min_lat = min([v.y for v in converted_vertices])
    max_lat = max([v.y for v in converted_vertices])

    # Convert boundary vertices to polygon for containment check as we're dealing with a hexagon (some grid points will be out of bounds)
    # This is also done to compute test data points within the hexagon.

    polygon = Polygon(converted_vertices)
    poly_df = gpd.GeoDataFrame(geometry=[polygon], crs=epsg)

    lons = np.arange(min_lon, max_lon, grid_spacing)
    lats = np.arange(min_lat, max_lat, grid_spacing)
    grid_points = np.array([(lon, lat) for lon in lons for lat in lats])

    valid_grid_points = [point for point in grid_points if Point(point).within(polygon)]

    grid_hex_x = np.array([point[0] for point in valid_grid_points])
    grid_hex_y = np.array([point[1] for point in valid_grid_points])
    
    return grid_hex_x, grid_hex_y

def generate_grid_over_data(data, epsg, grid_spacing, shapefile):
    """ Generate grid points over a region using its bounding box

    Args:
        data (gpd.GeoDataFrame): GeoDataFrame containing the data points (not needed anymore. using shapefile instead)
        epsg (int): EPSG code for the coordinate reference system
        grid_spacing (float): Spacing between grid points

    Returns:
        np.array: x-coordinates of grid points
        np.array: y-coordinates of grid points
    """
    bdry = gpd.read_file(shapefile, crs="EPSG:4326")
    bdry = bdry.to_crs(epsg)
    # Get bounding box
    min_lon, min_lat, max_lon, max_lat = bdry.total_bounds

    lons = np.arange(min_lon, max_lon, grid_spacing)
    lats = np.arange(min_lat, max_lat, grid_spacing)
    grid_points = np.array([(lon, lat) for lon in lons for lat in lats])

    valid_grid_points = [point for point in grid_points if Point(point).within(bdry.geometry[0])]

    grid_x = np.array([point[0] for point in valid_grid_points])
    grid_y = np.array([point[1] for point in valid_grid_points])
    
    return grid_x, grid_y