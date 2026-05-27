"""Evaluate clustering stability across sub-periods using pairwise ARI.

For each pair of sub-period cluster GeoJSONs, computes Adjusted Rand Index
over the intersection of hex cells present in both periods. Reports mean ARI
as a single stability score.

Usage:
    uv run python src/evaluate.py \\
        --city chicago \\
        --start 2024-01-01 --end 2024-12-31 \\
        --granularity week \\
        [--method idw] [--resolution 8]
"""

import json
import os
import sys
import glob
import click
import itertools
import numpy as np
import geopandas as gpd
from sklearn.metrics import adjusted_rand_score

SCRIBE_ROOT = "/data/taveesh/scribe"
OUT_DIR = os.path.join(SCRIBE_ROOT, "output")

_GRAN_FREQ = {"day": "D", "week": "W", "month": "MS"}


def _cluster_files(city: str, start: str, end: str, method: str, resolution: int) -> list[str]:
    """Return sorted list of cluster GeoJSON paths matching the city/method/resolution."""
    pattern = os.path.join(
        OUT_DIR,
        f"{city}_*_{method}_res{resolution}_clusters.geojson",
    )
    all_files = sorted(glob.glob(pattern))
    # Filter to those whose period falls within [start, end].
    import pandas as pd
    result = []
    for path in all_files:
        basename = os.path.basename(path)
        # Pattern: {city}_{period_start}_{period_end}_{method}_res{N}_clusters.geojson
        parts = basename.replace(f"{city}_", "").replace(f"_{method}_res{resolution}_clusters.geojson", "")
        # parts is now "{period_start}_{period_end}"
        try:
            ps, pe = parts.split("_", 1)
            if ps >= start and pe <= end:
                result.append(path)
        except ValueError:
            continue
    return result


def _load_labels(path: str) -> dict[str, int]:
    """Return {hex_id: cluster_label} for a cluster GeoJSON."""
    gdf = gpd.read_file(path)
    return dict(zip(gdf["hex_id"], gdf["cluster"].astype(int)))


def _pairwise_ari(
    labels_a: dict[str, int],
    labels_b: dict[str, int],
    min_overlap_frac: float = 0.5,
) -> float | None:
    """Compute ARI on the intersection of hex cells.

    Returns None if the intersection is smaller than min_overlap_frac of
    either period (unreliable comparison).
    """
    shared = set(labels_a) & set(labels_b)
    if len(shared) < min_overlap_frac * len(labels_a):
        return None
    if len(shared) < min_overlap_frac * len(labels_b):
        return None
    ya = [labels_a[h] for h in shared]
    yb = [labels_b[h] for h in shared]
    return adjusted_rand_score(ya, yb)


@click.command()
@click.option("--city", required=True)
@click.option("--start", required=True)
@click.option("--end", required=True)
@click.option("--granularity", required=True, type=click.Choice(["day", "week", "month"]))
@click.option("--method", default="idw", show_default=True, type=click.Choice(["idw", "loess", "kde"]))
@click.option("--resolution", default=8, show_default=True)
@click.option(
    "--min-overlap",
    default=0.5,
    show_default=True,
    help="Minimum fraction of hex cells that must overlap for a pair to be included.",
)
def main(city, start, end, granularity, method, resolution, min_overlap):
    files = _cluster_files(city, start, end, method, resolution)

    if len(files) < 2:
        click.echo(
            f"Found {len(files)} cluster file(s) — need at least 2 to compute ARI. "
            "Run 'make cluster' first.",
            err=True,
        )
        sys.exit(1)

    click.echo(f"Computing pairwise ARI across {len(files)} periods ...")

    all_labels = [_load_labels(f) for f in files]
    ari_values = []
    skipped = 0

    for (i, la), (j, lb) in itertools.combinations(enumerate(all_labels), 2):
        ari = _pairwise_ari(la, lb, min_overlap_frac=min_overlap)
        if ari is None:
            skipped += 1
        else:
            ari_values.append(ari)

    total_pairs = len(list(itertools.combinations(range(len(files)), 2)))

    if not ari_values:
        click.echo("No valid pairs (all had insufficient overlap). Try lowering --min-overlap.", err=True)
        sys.exit(1)

    result = {
        "city": city,
        "start": start,
        "end": end,
        "granularity": granularity,
        "method": method,
        "resolution": resolution,
        "n_periods": len(files),
        "n_pairs_total": total_pairs,
        "n_pairs_used": len(ari_values),
        "n_pairs_skipped": skipped,
        "stability_mean_ari": float(np.mean(ari_values)),
        "stability_std_ari": float(np.std(ari_values)),
        "stability_min_ari": float(np.min(ari_values)),
        "stability_max_ari": float(np.max(ari_values)),
    }

    out_path = os.path.join(
        OUT_DIR, f"{city}_{start}_{end}_{granularity}_{method}_stability.json"
    )
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    click.echo(
        f"\nStability (mean ARI): {result['stability_mean_ari']:.4f} "
        f"± {result['stability_std_ari']:.4f}  "
        f"[min={result['stability_min_ari']:.4f}, max={result['stability_max_ari']:.4f}]"
    )
    click.echo(f"  pairs used: {len(ari_values)}/{total_pairs}  (skipped {skipped} low-overlap pairs)")
    click.echo(f"  saved → {out_path}")


if __name__ == "__main__":
    main()
