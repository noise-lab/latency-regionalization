"""Fetch M-Lab NDT latency measurements from BigQuery.

Caches results locally as Parquet so re-running never re-queries BigQuery.

Usage:
    uv run python src/fetch_data.py \\
        --city chicago \\
        --start 2024-01-01 --end 2024-12-31 \\
        [--dry-run] [--accuracy-radius 50]

Environment variables:
    GOOGLE_APPLICATION_CREDENTIALS  Path to service account JSON key file.
    SAMPLE_PCT                       If set, adds TABLESAMPLE SYSTEM(N PERCENT)
                                     for development use.
"""

import os
import sys
import click
import pandas as pd
import geopandas as gpd
from google.cloud import bigquery
from google.oauth2 import service_account

from cities import get_city_bbox, get_city_boundary
from data_utils import save_parquet

SCRIBE_ROOT = "/data/taveesh/scribe"
RAW_DIR = os.path.join(SCRIBE_ROOT, "raw")

_QUERY_TEMPLATE = """
SELECT
  a.TestTime                        AS test_time,
  a.MinRTT                          AS latency_ms,
  client.Geo.Latitude               AS latitude,
  client.Geo.Longitude              AS longitude,
  client.Geo.AccuracyRadius         AS accuracy_radius,
  client.Network.ASName             AS attr_provider_name,
  client.Network.ASNumber           AS asn
FROM `measurement-lab.ndt.unified_downloads`{sample_clause}
WHERE DATE(a.TestTime) BETWEEN @start_date AND @end_date
  AND client.Geo.Latitude  BETWEEN @min_lat  AND @max_lat
  AND client.Geo.Longitude BETWEEN @min_lon  AND @max_lon
  AND client.Geo.AccuracyRadius <= @accuracy_radius
  AND a.MinRTT IS NOT NULL
  AND a.MinRTT > 0
"""


def _cache_path(city: str, start: str, end: str) -> str:
    return os.path.join(RAW_DIR, f"{city}_{start}_{end}.parquet")


def _build_client() -> bigquery.Client:
    key_path = os.environ.get("GOOGLE_APPLICATION_CREDENTIALS")
    if not key_path:
        raise EnvironmentError(
            "GOOGLE_APPLICATION_CREDENTIALS is not set. "
            "Point it to your service account JSON key file."
        )
    credentials = service_account.Credentials.from_service_account_file(
        key_path,
        scopes=["https://www.googleapis.com/auth/bigquery"],
    )
    return bigquery.Client(credentials=credentials, project=credentials.project_id)


def _run_query(
    city: str,
    start: str,
    end: str,
    bbox: dict,
    accuracy_radius: int,
    dry_run: bool,
    sample_pct: float | None,
) -> pd.DataFrame | None:
    sample_clause = f"\nTABLESAMPLE SYSTEM ({sample_pct} PERCENT)" if sample_pct else ""
    query = _QUERY_TEMPLATE.format(sample_clause=sample_clause)

    job_config = bigquery.QueryJobConfig(
        dry_run=dry_run,
        use_query_cache=True,
        query_parameters=[
            bigquery.ScalarQueryParameter("start_date", "DATE", start),
            bigquery.ScalarQueryParameter("end_date", "DATE", end),
            bigquery.ScalarQueryParameter("min_lat", "FLOAT64", bbox["min_lat"]),
            bigquery.ScalarQueryParameter("max_lat", "FLOAT64", bbox["max_lat"]),
            bigquery.ScalarQueryParameter("min_lon", "FLOAT64", bbox["min_lon"]),
            bigquery.ScalarQueryParameter("max_lon", "FLOAT64", bbox["max_lon"]),
            bigquery.ScalarQueryParameter("accuracy_radius", "INT64", accuracy_radius),
        ],
    )

    client = _build_client()
    job = client.query(query, job_config=job_config)

    if dry_run:
        gb = job.total_bytes_processed / 1e9
        print(f"[dry-run] Estimated bytes: {gb:.2f} GB")
        return None

    print(f"Running BigQuery job for {city} ({start} → {end}) ...")
    df = job.result().to_dataframe(progress_bar_type="tqdm")
    print(f"  fetched {len(df):,} rows.")
    return df


def _clip_to_city(df: pd.DataFrame, city: str) -> pd.DataFrame:
    """Drop measurements that fall outside the city boundary polygon."""
    boundary = get_city_boundary(city)
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["longitude"], df["latitude"]),
        crs="EPSG:4326",
    )
    polygon = boundary.geometry.iloc[0]
    mask = gdf.geometry.within(polygon)
    kept = int(mask.sum())
    dropped = len(df) - kept
    if dropped:
        print(f"  containment clip: dropped {dropped:,} rows outside {city} boundary.")
    return df[mask].drop(columns="geometry", errors="ignore").reset_index(drop=True)


@click.command()
@click.option("--city", required=True, help="City name (must exist in data/cities.geojson).")
@click.option("--start", required=True, help="Start date YYYY-MM-DD.")
@click.option("--end", required=True, help="End date YYYY-MM-DD.")
@click.option("--accuracy-radius", default=50, show_default=True, help="Max geolocation accuracy radius (km).")
@click.option("--dry-run", is_flag=True, default=False, help="Estimate bytes without running the query.")
def main(city: str, start: str, end: str, accuracy_radius: int, dry_run: bool) -> None:
    out_path = _cache_path(city, start, end)

    if not dry_run and os.path.exists(out_path):
        print(f"Cache hit: {out_path} — skipping BigQuery.")
        return

    bbox = get_city_bbox(city)
    sample_pct_env = os.environ.get("SAMPLE_PCT")
    sample_pct = float(sample_pct_env) if sample_pct_env else None

    df = _run_query(city, start, end, bbox, accuracy_radius, dry_run, sample_pct)

    if dry_run or df is None:
        return

    df = _clip_to_city(df, city)

    save_parquet(df, out_path)
    print(f"Saved {len(df):,} rows → {out_path}")


if __name__ == "__main__":
    main()
