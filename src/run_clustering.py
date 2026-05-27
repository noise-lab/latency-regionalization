"""Run SKATER spatial clustering on hex-aggregated latency data.

Usage:
    uv run python src/run_clustering.py \\
        --city chicago \\
        --period-start 2024-01-01 --period-end 2024-01-07 \\
        [--method idw] [--resolution 8] \\
        [--n-clusters auto] [--distance-metric Euclidean] [--n-bootstraps 100]
"""

import os
import sys
import click
import geopandas as gpd

from clustering import SKATERClustering

SCRIBE_ROOT = "/data/taveesh/scribe"
AGG_DIR = os.path.join(SCRIBE_ROOT, "aggregated")
OUT_DIR = os.path.join(SCRIBE_ROOT, "output")

# Aggregate columns used as features for SKATER clustering.
CLUSTER_FEATURES = ["mean", "p10", "p50", "p90", "latency_reduction", "inequality_ratio"]

DISTANCE_METRICS = ["Manhattan", "Euclidean", "Chebyshev", "Cosine", "Correlation", "Canberra"]


def _agg_path(city, period_start, period_end, method, resolution):
    return os.path.join(AGG_DIR, f"{city}_{period_start}_{period_end}_{method}_res{resolution}.parquet")


def _out_path(city, period_start, period_end, method, resolution):
    return os.path.join(OUT_DIR, f"{city}_{period_start}_{period_end}_{method}_res{resolution}_clusters.geojson")


@click.command()
@click.option("--city", required=True)
@click.option("--period-start", required=True)
@click.option("--period-end", required=True)
@click.option("--method", default="idw", show_default=True, type=click.Choice(["idw", "loess", "kde"]))
@click.option("--resolution", default=8, show_default=True)
@click.option(
    "--n-clusters",
    default="auto",
    show_default=True,
    help="Number of clusters, or 'auto' to detect via get_optimal_clusters().",
)
@click.option(
    "--distance-metric",
    default="Euclidean",
    show_default=True,
    type=click.Choice(DISTANCE_METRICS),
)
@click.option("--n-bootstraps", default=100, show_default=True)
def main(city, period_start, period_end, method, resolution, n_clusters, distance_metric, n_bootstraps):
    agg_path = _agg_path(city, period_start, period_end, method, resolution)
    if not os.path.exists(agg_path):
        click.echo(f"Aggregated data not found: {agg_path}\nRun 'make aggregate' first.", err=True)
        sys.exit(1)

    out_path = _out_path(city, period_start, period_end, method, resolution)
    if os.path.exists(out_path):
        click.echo(f"Cache hit: {out_path} — skipping.")
        return

    gdf = gpd.read_parquet(agg_path)

    # Filter to features that exist in the GeoDataFrame.
    available_features = [f for f in CLUSTER_FEATURES if f in gdf.columns]
    if not available_features:
        click.echo("No cluster features found in aggregated data.", err=True)
        sys.exit(1)

    click.echo(
        f"[{city} {period_start}→{period_end}] SKATER on {len(gdf)} hex cells "
        f"({len(available_features)} features, metric={distance_metric})"
    )

    # Input is already hex-aggregated; set use_hexagons=False to skip re-tessellation.
    sc = SKATERClustering(
        data=gdf,
        resolution=resolution,
        aggregates=available_features,
        metric="mean",
        use_hexagons=False,
        n_bootstraps=n_bootstraps,
    )

    if n_clusters == "auto":
        click.echo("  detecting optimal cluster count ...")
        k, floor = sc.get_optimal_clusters(distance_metric)
        click.echo(f"  optimal K={k}, floor={floor}")
    else:
        k = int(n_clusters)
        num_units = len(gdf)
        floor = num_units // 10

    result = sc.cluster(k, floor, distance_metric)

    os.makedirs(OUT_DIR, exist_ok=True)
    result.to_file(out_path, driver="GeoJSON")
    click.echo(f"  saved {k} clusters → {out_path}")


if __name__ == "__main__":
    main()
