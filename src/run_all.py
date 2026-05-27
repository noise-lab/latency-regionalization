"""Orchestrator: run a pipeline stage for every sub-period in a date range.

Called by Makefile targets (interpolate, aggregate, cluster). Splits the
date range into sub-periods according to GRANULARITY, then invokes the
corresponding script for each period.

Usage (via Makefile, not directly):
    uv run python src/run_all.py \\
        --city chicago --start 2024-01-01 --end 2024-12-31 \\
        --granularity week --stage interpolate \\
        --method idw --resolution 8 ...
"""

import subprocess
import sys
import os
import click
import pandas as pd

_FREQ = {"day": "D", "week": "W-MON", "month": "MS"}

_STAGE_SCRIPT = {
    "interpolate": "run_interpolation.py",
    "aggregate": "aggregate.py",
    "cluster": "run_clustering.py",
}

SRC_DIR = os.path.dirname(__file__)


def _sub_periods(start: str, end: str, granularity: str) -> list[tuple[str, str]]:
    """Return list of (period_start, period_end) strings covering [start, end]."""
    freq = _FREQ[granularity]
    starts = pd.date_range(start=start, end=end, freq=freq)
    periods = []
    for i, ps in enumerate(starts):
        pe = starts[i + 1] - pd.Timedelta(days=1) if i + 1 < len(starts) else pd.Timestamp(end)
        pe = min(pe, pd.Timestamp(end))
        ps_str = ps.strftime("%Y-%m-%d")
        pe_str = pe.strftime("%Y-%m-%d")
        if ps_str <= end:
            periods.append((ps_str, pe_str))
    # Edge case: if start is between two period boundaries, include a leading partial period.
    if not periods or periods[0][0] > start:
        first_end = (pd.Timestamp(periods[0][0]) - pd.Timedelta(days=1)).strftime("%Y-%m-%d") if periods else end
        periods.insert(0, (start, first_end))
    return periods


def _run_stage(script: str, common_args: list[str], extra_args: list[str]) -> int:
    cmd = ["uv", "run", "python", os.path.join(SRC_DIR, script)] + common_args + extra_args
    result = subprocess.run(cmd)
    return result.returncode


@click.command(context_settings={"ignore_unknown_options": True, "allow_extra_args": True})
@click.option("--city", required=True)
@click.option("--start", required=True)
@click.option("--end", required=True)
@click.option("--granularity", required=True, type=click.Choice(["day", "week", "month"]))
@click.option("--stage", required=True, type=click.Choice(["interpolate", "aggregate", "cluster"]))
@click.option("--method", default="idw")
@click.option("--resolution", default=8)
@click.option("--min-samples", default=5)
@click.option("--n-clusters", default="auto")
@click.option("--distance-metric", default="Euclidean")
@click.option("--n-bootstraps", default=100)
# Interpolation-specific options forwarded as-is:
@click.option("--grid-spacing", default=500.0)
@click.option("--p", default=2.0)
@click.option("--span", default=0.3)
@click.option("--degree", default=2)
@click.option("--c", default=0.5)
@click.option("--k", default=10)
def main(
    city, start, end, granularity, stage,
    method, resolution, min_samples, n_clusters, distance_metric, n_bootstraps,
    grid_spacing, p, span, degree, c, k,
):
    periods = _sub_periods(start, end, granularity)
    click.echo(
        f"[{city}] {stage} | granularity={granularity} | "
        f"{len(periods)} periods | method={method} | res={resolution}"
    )

    script = _STAGE_SCRIPT[stage]
    failed = 0

    for period_start, period_end in periods:
        click.echo(f"  → {period_start} to {period_end}")

        common = [
            "--city", city,
            "--period-start", period_start,
            "--period-end", period_end,
            "--method", method,
            "--resolution", str(resolution),
        ]

        if stage == "interpolate":
            extra = [
                "--full-start", start,
                "--full-end", end,
                "--grid-spacing", str(grid_spacing),
                "--p", str(p),
                "--span", str(span),
                "--degree", str(degree),
                "--c", str(c),
                "--k", str(k),
            ]
        elif stage == "aggregate":
            extra = ["--min-samples", str(min_samples)]
        elif stage == "cluster":
            extra = [
                "--n-clusters", str(n_clusters),
                "--distance-metric", distance_metric,
                "--n-bootstraps", str(n_bootstraps),
            ]
        else:
            extra = []

        rc = _run_stage(script, common, extra)
        if rc != 0:
            click.echo(f"    FAILED (exit {rc})", err=True)
            failed += 1

    if failed:
        click.echo(f"\n{failed}/{len(periods)} periods failed.", err=True)
        sys.exit(1)

    click.echo(f"\nCompleted {len(periods)} periods.")


if __name__ == "__main__":
    main()
