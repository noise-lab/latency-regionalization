"""City boundary and metadata lookups backed by data/cities.geojson.

Run src/seed_cities.py once to populate the GeoJSON before using these
functions in the pipeline.
"""

import json
import os
import geopandas as gpd
from shapely.geometry import shape

_CITIES_FILE = os.path.join(os.path.dirname(__file__), "..", "data", "cities.geojson")

# Approximate UTM EPSG codes for major US regions. The lookup is coarse — for
# high-precision work a proper UTM zone calculation should replace this table.
_APPROX_UTM = {
    # lon range          EPSG
    (-180, -162): 32601,
    (-162, -156): 32602,
    (-156, -150): 32603,
    (-150, -144): 32604,
    (-144, -138): 32605,
    (-138, -132): 32606,
    (-132, -126): 32607,
    (-126, -120): 32608,
    (-120, -114): 32609,
    (-114, -108): 32610,
    (-108, -102): 32611,
    (-102,  -96): 32612,
    ( -96,  -90): 32613,
    ( -90,  -84): 32614,
    ( -84,  -78): 32615,
    ( -78,  -72): 32616,
    ( -72,  -66): 32617,
    ( -66,  -60): 32618,
}


def _load_collection() -> list[dict]:
    with open(_CITIES_FILE) as f:
        return json.load(f)["features"]


def _find_feature(city_name: str) -> dict:
    key = city_name.strip().lower().replace(" ", "_")
    for feat in _load_collection():
        if feat["properties"].get("city_key") == key:
            return feat
    raise KeyError(
        f"City '{city_name}' not found in data/cities.geojson. "
        "Run 'uv run python src/seed_cities.py' to populate it."
    )


def get_city_bbox(city_name: str) -> dict:
    """Return the bounding box for a city.

    Returns
    -------
    dict with keys: min_lat, max_lat, min_lon, max_lon
    """
    feat = _find_feature(city_name)
    props = feat["properties"]
    return {
        "min_lat": props["min_lat"],
        "max_lat": props["max_lat"],
        "min_lon": props["min_lon"],
        "max_lon": props["max_lon"],
    }


def get_city_boundary(city_name: str) -> gpd.GeoDataFrame:
    """Return a single-row GeoDataFrame with the city boundary polygon."""
    feat = _find_feature(city_name)
    geometry = shape(feat["geometry"])
    return gpd.GeoDataFrame(
        [feat["properties"]],
        geometry=[geometry],
        crs="EPSG:4326",
    )


def get_city_utm_epsg(city_name: str) -> int:
    """Return an appropriate projected (UTM) EPSG code for the city."""
    bbox = get_city_bbox(city_name)
    center_lon = (bbox["min_lon"] + bbox["max_lon"]) / 2
    for (lo, hi), epsg in _APPROX_UTM.items():
        if lo <= center_lon < hi:
            return epsg
    # Fall back to Web Mercator for non-US cities or edge cases.
    return 3857
