"""Overlay H3 hexagons on an interpolated grid and compute per-hex aggregates.

Takes interpolated grid Parquet (grid_x, grid_y, z_smooth) and produces a
GeoParquet with one row per hex cell containing latency distribution stats.

Usage:
    uv run python src/aggregate.py \\
        --city chicago \\
        --period-start 2024-01-01 --period-end 2024-01-07 \\
        [--method idw] [--resolution 8] [--min-samples 5]
"""

import os
import sys
import click
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

from cities import get_city_utm_epsg
from data_utils import load_parquet
from geo_utils import assign_hex_ids
from clustering import SpatialClustering

SCRIBE_ROOT = "/data/taveesh/scribe"
INTERP_DIR = os.path.join(SCRIBE_ROOT, "interpolated")
AGG_DIR = os.path.join(SCRIBE_ROOT, "aggregated")


def _interp_path(city: str, period_start: str, period_end: str, method: str) -> str:
    return os.path.join(INTERP_DIR, f"{city}_{period_start}_{period_end}_{method}.parquet")


def _out_path(city: str, period_start: str, period_end: str, method: str, resolution: int) -> str:
    return os.path.join(AGG_DIR, f"{city}_{period_start}_{period_end}_{method}_res{resolution}.parquet")


@click.command()
@click.option("--city", required=True)
@click.option("--period-start", required=True)
@click.option("--period-end", required=True)
@click.option("--method", default="idw", show_default=True, type=click.Choice(["idw", "loess", "kde"]))
@click.option("--resolution", default=8, show_default=True, help="H3 resolution for hexagon tessellation.")
@click.option("--min-samples", default=5, show_default=True, help="Minimum grid points per hex cell.")
def main(city, period_start, period_end, method, resolution, min_samples):
    interp_path = _interp_path(city, period_start, period_end, method)
    if not os.path.exists(interp_path):
        click.echo(f"Interpolated data not found: {interp_path}\nRun 'make interpolate' first.", err=True)
        sys.exit(1)

    out_path = _out_path(city, period_start, period_end, method, resolution)
    if os.path.exists(out_path):
        click.echo(f"Cache hit: {out_path} — skipping.")
        return

    df = load_parquet(interp_path)
    df = df.dropna(subset=["z_smooth"])

    if df.empty:
        click.echo(f"No valid interpolated data in {interp_path}. Skipping.", err=True)
        return

    epsg = get_city_utm_epsg(city)

    # Build GeoDataFrame in the projected CRS used during interpolation, then
    # reproject to EPSG:4326 for H3 assignment (H3 requires lat/lon).
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["grid_x"], df["grid_y"]),
        crs=epsg,
    ).to_crs("EPSG:4326")

    # Rename the interpolated value to a column we can aggregate on.
    gdf["latency_ms"] = gdf["z_smooth"]

    click.echo(f"[{city} {period_start}→{period_end}] assigning H3 res-{resolution} cells ...")
    gdf, tessellation = assign_hex_ids(gdf, resolution)

    # Drop sparse cells.
    cell_counts = gdf.groupby("hex_id").size()
    valid_cells = cell_counts[cell_counts >= min_samples].index
    dropped = (~gdf["hex_id"].isin(valid_cells)).sum()
    if dropped:
        click.echo(f"  dropped {dropped:,} grid points in cells with <{min_samples} samples.")
    gdf = gdf[gdf["hex_id"].isin(valid_cells)].reset_index(drop=True)
    tessellation = tessellation[tessellation["hex_id"].isin(valid_cells)].reset_index(drop=True)

    if gdf.empty:
        click.echo("  no cells survived the min-samples filter. Skipping.", err=True)
        return

    # Reuse SpatialClustering.aggregate() to compute per-hex distribution stats.
    dummy_clustering = SpatialClustering(
        data=gdf,
        resolution=resolution,
        aggregates=["mean"],
        metric="latency_ms",
        use_hexagons=False,
    )
    hex_agg = dummy_clustering.aggregate(gdf, over="hex_id", tessellation=tessellation)

    os.makedirs(AGG_DIR, exist_ok=True)
    hex_agg.to_parquet(out_path)
    click.echo(f"  saved {len(hex_agg):,} hex cells → {out_path}")


if __name__ == "__main__":
    main()
