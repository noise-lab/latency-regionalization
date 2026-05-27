"""One-time script to populate data/cities.geojson from OSM Nominatim.

Usage:
    uv run python src/seed_cities.py

Adds entries for each city name listed in CITIES. Safe to re-run — existing
entries are skipped and only missing cities are fetched.
"""

import json
import os
import time
import requests

CITIES_FILE = os.path.join(os.path.dirname(__file__), "..", "data", "cities.geojson")

# Extend this list to add support for more cities.
CITIES = [
    "Chicago, Illinois, USA",
    "New York City, New York, USA",
    "Los Angeles, California, USA",
    "Houston, Texas, USA",
    "Phoenix, Arizona, USA",
    "Philadelphia, Pennsylvania, USA",
    "San Antonio, Texas, USA",
    "San Diego, California, USA",
    "Dallas, Texas, USA",
    "San Jose, California, USA",
    "Austin, Texas, USA",
    "Jacksonville, Florida, USA",
    "Fort Worth, Texas, USA",
    "Columbus, Ohio, USA",
    "Charlotte, North Carolina, USA",
    "Indianapolis, Indiana, USA",
    "San Francisco, California, USA",
    "Seattle, Washington, USA",
    "Denver, Colorado, USA",
    "Nashville, Tennessee, USA",
    "Oklahoma City, Oklahoma, USA",
    "El Paso, Texas, USA",
    "Washington, District of Columbia, USA",
    "Las Vegas, Nevada, USA",
    "Louisville, Kentucky, USA",
    "Memphis, Tennessee, USA",
    "Portland, Oregon, USA",
    "Baltimore, Maryland, USA",
    "Milwaukee, Wisconsin, USA",
    "Albuquerque, New Mexico, USA",
    "Tucson, Arizona, USA",
    "Fresno, California, USA",
    "Mesa, Arizona, USA",
    "Kansas City, Missouri, USA",
    "Atlanta, Georgia, USA",
    "Omaha, Nebraska, USA",
    "Colorado Springs, Colorado, USA",
    "Raleigh, North Carolina, USA",
    "Long Beach, California, USA",
    "Virginia Beach, Virginia, USA",
    "Minneapolis, Minnesota, USA",
    "Tampa, Florida, USA",
    "New Orleans, Louisiana, USA",
    "Arlington, Texas, USA",
    "Wichita, Kansas, USA",
    "Bakersfield, California, USA",
    "Aurora, Colorado, USA",
    "Anaheim, California, USA",
    "Santa Ana, California, USA",
    "Corpus Christi, Texas, USA",
    "Riverside, California, USA",
    "Lexington, Kentucky, USA",
    "Stockton, California, USA",
    "Pittsburgh, Pennsylvania, USA",
    "St. Paul, Minnesota, USA",
    "Anchorage, Alaska, USA",
    "Cincinnati, Ohio, USA",
    "Greensboro, North Carolina, USA",
    "Toledo, Ohio, USA",
    "Newark, New Jersey, USA",
    "Plano, Texas, USA",
    "Henderson, Nevada, USA",
    "Orlando, Florida, USA",
    "St. Louis, Missouri, USA",
    "Laredo, Texas, USA",
    "Norfolk, Virginia, USA",
    "Madison, Wisconsin, USA",
    "Durham, North Carolina, USA",
    "Lubbock, Texas, USA",
    "Winston-Salem, North Carolina, USA",
    "Garland, Texas, USA",
    "Glendale, Arizona, USA",
    "Hialeah, Florida, USA",
    "Reno, Nevada, USA",
    "Baton Rouge, Louisiana, USA",
    "Irvine, California, USA",
    "Chesapeake, Virginia, USA",
    "Irving, Texas, USA",
    "Scottsdale, Arizona, USA",
    "North Las Vegas, Nevada, USA",
    "Fremont, California, USA",
    "Gilbert, Arizona, USA",
    "San Bernardino, California, USA",
    "Birmingham, Alabama, USA",
    "Boise, Idaho, USA",
    "Rochester, New York, USA",
]

NOMINATIM_URL = "https://nominatim.openstreetmap.org/search"
HEADERS = {"User-Agent": "latency-regionalization/1.0 (research project)"}


def _canonical_name(city_query: str) -> str:
    """Extract a short canonical key from a full query string."""
    return city_query.split(",")[0].strip().lower().replace(" ", "_")


def fetch_boundary(city_query: str) -> dict | None:
    """Fetch the city polygon from Nominatim and return a GeoJSON feature."""
    params = {
        "q": city_query,
        "format": "geojson",
        "polygon_geojson": 1,
        "limit": 1,
    }
    resp = requests.get(NOMINATIM_URL, params=params, headers=HEADERS, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    if not data.get("features"):
        return None
    feature = data["features"][0]
    feature["properties"]["city_key"] = _canonical_name(city_query)
    feature["properties"]["query"] = city_query
    bbox = feature.get("bbox")
    if bbox:
        feature["properties"]["min_lon"] = bbox[0]
        feature["properties"]["min_lat"] = bbox[1]
        feature["properties"]["max_lon"] = bbox[2]
        feature["properties"]["max_lat"] = bbox[3]
    return feature


def main() -> None:
    with open(CITIES_FILE) as f:
        collection = json.load(f)

    existing_keys = {
        feat["properties"].get("city_key")
        for feat in collection.get("features", [])
    }

    added = 0
    for city_query in CITIES:
        key = _canonical_name(city_query)
        if key in existing_keys:
            print(f"  skip (exists): {key}")
            continue
        print(f"  fetching: {city_query} ...", end=" ", flush=True)
        try:
            feature = fetch_boundary(city_query)
            if feature is None:
                print("NOT FOUND")
            else:
                collection["features"].append(feature)
                existing_keys.add(key)
                added += 1
                print("ok")
        except Exception as e:
            print(f"ERROR: {e}")
        time.sleep(1)  # Nominatim rate limit: 1 request/second

    with open(CITIES_FILE, "w") as f:
        json.dump(collection, f, separators=(",", ":"))

    print(f"\nAdded {added} cities. Total: {len(collection['features'])}.")


if __name__ == "__main__":
    main()
