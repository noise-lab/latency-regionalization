"""Interpolate raw measurements to a regular grid for a single time period.

Usage:
    uv run python src/run_interpolation.py \\
        --city chicago \\
        --period-start 2024-01-01 --period-end 2024-01-07 \\
        [--method idw] [--resolution 8] [--p 2] [--span 0.3] [--degree 2] \\
        [--c 0.5] [--k 10] [--grid-spacing 500]
"""

import os
import sys
import click
import numpy as np
import pandas as pd

from cities import get_city_boundary, get_city_utm_epsg
from data_utils import load_parquet, save_parquet
from interpolation import IDW, LOESS_Interpolator, STBKR_Interpolator

SCRIBE_ROOT = "/data/taveesh/scribe"
RAW_DIR = os.path.join(SCRIBE_ROOT, "raw")
INTERP_DIR = os.path.join(SCRIBE_ROOT, "interpolated")


def _raw_path(city: str, start: str, end: str) -> str:
    return os.path.join(RAW_DIR, f"{city}_{start}_{end}.parquet")


def _out_path(city: str, period_start: str, period_end: str, method: str) -> str:
    return os.path.join(INTERP_DIR, f"{city}_{period_start}_{period_end}_{method}.parquet")


def _build_interpolator(method: str, data: pd.DataFrame, params: dict, shapefile: str, **kwargs):
    if method == "idw":
        return IDW(data, params, p=kwargs["p"], shapefile=shapefile)
    elif method == "loess":
        return LOESS_Interpolator(
            data, params,
            span=kwargs["span"],
            degree=kwargs["degree"],
            shapefile=shapefile,
        )
    elif method == "kde":
        return STBKR_Interpolator(
            data, params,
            c=kwargs["c"],
            k=kwargs["k"],
            shapefile=shapefile,
        )
    else:
        raise ValueError(f"Unknown interpolation method: {method!r}. Choose idw, loess, or kde.")


@click.command()
@click.option("--city", required=True)
@click.option("--period-start", required=True, help="Period start YYYY-MM-DD.")
@click.option("--period-end", required=True, help="Period end YYYY-MM-DD.")
@click.option("--full-start", default=None, help="Full data range start (defaults to period-start).")
@click.option("--full-end", default=None, help="Full data range end (defaults to period-end).")
@click.option("--method", default="idw", show_default=True, type=click.Choice(["idw", "loess", "kde"]))
@click.option("--resolution", default=8, show_default=True, help="H3 resolution (unused in region-wide mode but stored in params).")
@click.option("--grid-spacing", default=500.0, show_default=True, help="Grid spacing in projected CRS units (metres).")
@click.option("--p", default=2.0, show_default=True, help="IDW power parameter.")
@click.option("--span", default=0.3, show_default=True, help="LOESS span (fraction of data).")
@click.option("--degree", default=2, show_default=True, help="LOESS polynomial degree.")
@click.option("--c", default=0.5, show_default=True, help="KDE bandwidth scaling constant.")
@click.option("--k", default=10, show_default=True, help="KDE number of nearest neighbours.")
def main(
    city, period_start, period_end, full_start, full_end,
    method, resolution, grid_spacing, p, span, degree, c, k,
):
    full_start = full_start or period_start
    full_end = full_end or period_end

    raw_path = _raw_path(city, full_start, full_end)
    if not os.path.exists(raw_path):
        click.echo(f"Raw data not found: {raw_path}\nRun 'make fetch' first.", err=True)
        sys.exit(1)

    out_path = _out_path(city, period_start, period_end, method)
    if os.path.exists(out_path):
        click.echo(f"Cache hit: {out_path} — skipping.")
        return

    df = load_parquet(raw_path)
    # Slice to the requested sub-period.
    df["test_time"] = pd.to_datetime(df["test_time"])
    mask = (df["test_time"].dt.date >= pd.Timestamp(period_start).date()) & \
           (df["test_time"].dt.date <= pd.Timestamp(period_end).date())
    df = df[mask].reset_index(drop=True)

    if df.empty:
        click.echo(f"No data in period {period_start} → {period_end}. Skipping.", err=True)
        return

    click.echo(f"[{city} {period_start}→{period_end}] {method}: {len(df):,} measurements")

    epsg = get_city_utm_epsg(city)
    boundary = get_city_boundary(city)
    # Write boundary to a temp file so the interpolator can read it as a shapefile.
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        shapefile = os.path.join(tmpdir, "boundary.gpkg")
        boundary.to_file(shapefile, driver="GPKG")

        params = {
            "coord_columns": ["longitude", "latitude"],
            "h3_resolution": resolution,
            "agg_fns": {"mean": np.mean},
            "grid_spacing": grid_spacing,
            "metric": "latency_ms",
            "epsg": epsg,
        }

        interp = _build_interpolator(
            method, df, params, shapefile,
            p=p, span=span, degree=degree, c=c, k=k,
        )
        result = interp.fit()

    result = result.dropna(subset=["z_smooth"])
    save_parquet(result, out_path)
    click.echo(f"  saved {len(result):,} grid points → {out_path}")


if __name__ == "__main__":
    main()
