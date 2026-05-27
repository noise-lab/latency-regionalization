import h3
import geopandas as gpd
from shapely.geometry import Polygon, Point
import numpy as np


def _boundary_to_polygon(cell: str) -> Polygon:
    """Convert an H3 cell to a Shapely Polygon.

    h3.cell_to_boundary() returns (lat, lng) pairs; Shapely expects (lng, lat).
    """
    return Polygon([(lng, lat) for lat, lng in h3.cell_to_boundary(cell)])


def create_hexagon_around_point(h3_resolution: int, point: Point) -> Polygon:
    cell = h3.latlng_to_cell(point.y, point.x, h3_resolution)
    return _boundary_to_polygon(cell)


def assign_hex_ids(
    data: gpd.GeoDataFrame, h3_resolution: int
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Assign H3 cell IDs to each point in data.

    Returns:
        data: input GeoDataFrame with added hex_id and hex_centroid columns.
        hdf: GeoDataFrame of unique hexagon polygons keyed by hex_id.
    """
    cells = [
        h3.latlng_to_cell(row.geometry.y, row.geometry.x, h3_resolution)
        for _, row in data.iterrows()
    ]
    data = data.copy()
    data["hex_id"] = cells
    data["hex_centroid"] = data["hex_id"].apply(
        lambda c: Point(h3.cell_to_latlng(c)[1], h3.cell_to_latlng(c)[0])
    )
    hdf = gpd.GeoDataFrame(
        {"hex_id": cells},
        geometry=[_boundary_to_polygon(c) for c in cells],
        crs="EPSG:4326",
    )
    return data, hdf


def containment_filter(h3_cell: str, points: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Return points that fall within an H3 cell."""
    polygon = _boundary_to_polygon(h3_cell)
    points = points.to_crs("EPSG:4326")
    return points[points.geometry.apply(lambda p: p.within(polygon))]


def containment_filter_generic(
    bdry: gpd.GeoDataFrame, points: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Return points that fall within a polygon boundary GeoDataFrame."""
    polygon = bdry.geometry.iloc[0]
    points = points.to_crs(bdry.crs)
    return points[points.geometry.apply(lambda p: p.within(polygon))]


def generate_grid_points(
    h3_cell: str, epsg: int, grid_spacing: float
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a grid of points inside an H3 cell in a projected CRS."""
    boundary_vertices = [(lng, lat) for lat, lng in h3.cell_to_boundary(h3_cell)]
    bdf = gpd.GeoDataFrame(
        geometry=[Point(v) for v in boundary_vertices], crs="EPSG:4326"
    ).to_crs(epsg)
    converted_vertices = list(bdf.geometry)

    min_x = min(v.x for v in converted_vertices)
    max_x = max(v.x for v in converted_vertices)
    min_y = min(v.y for v in converted_vertices)
    max_y = max(v.y for v in converted_vertices)

    polygon = Polygon([(v.x, v.y) for v in converted_vertices])

    xs = np.arange(min_x, max_x, grid_spacing)
    ys = np.arange(min_y, max_y, grid_spacing)
    grid = [(x, y) for x in xs for y in ys if Point(x, y).within(polygon)]

    if not grid:
        return np.array([]), np.array([])

    grid_x = np.array([p[0] for p in grid])
    grid_y = np.array([p[1] for p in grid])
    return grid_x, grid_y


def generate_grid_over_data(
    epsg: int, grid_spacing: float, shapefile: str
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a grid of points over a region defined by a shapefile."""
    bdry = gpd.read_file(shapefile).to_crs(epsg)
    min_x, min_y, max_x, max_y = bdry.total_bounds
    polygon = bdry.geometry.iloc[0]

    xs = np.arange(min_x, max_x, grid_spacing)
    ys = np.arange(min_y, max_y, grid_spacing)
    grid = [(x, y) for x in xs for y in ys if Point(x, y).within(polygon)]

    if not grid:
        return np.array([]), np.array([])

    grid_x = np.array([p[0] for p in grid])
    grid_y = np.array([p[1] for p in grid])
    return grid_x, grid_y
