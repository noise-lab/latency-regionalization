# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**SCRIBE** (Spatially-Constrained Regionalization for Inference of Broadband Equity) — research code for "Beyond Data Points: Regionalizing Crowdsourced Latency Measurements" (ACM SIGMETRICS 2024). Implements an end-to-end pipeline that ingests M-Lab NDT latency data from BigQuery and produces spatial cluster maps per city and time period.

**Pipeline stages** (in order):
1. **Fetch** — query M-Lab BigQuery for `MinRTT` measurements within a city boundary and date range; cache as Parquet.
2. **Interpolate** — fit IDW / LOESS / KDE to raw point measurements and produce a regular interpolated grid.
3. **Aggregate** — overlay H3 hexagons on the grid; compute per-hex latency distribution stats.
4. **Cluster** — run SKATER spatial clustering on hex aggregates.
5. **Evaluate** — compute pairwise Adjusted Rand Index across sub-period clusters; report stability score.

All intermediates and outputs are stored under `/data/taveesh/scribe/`.

## Dependencies and Environment

Uses **uv** for all dependency management. Never use pip directly.

```bash
make setup        # runs uv sync; creates .venv and uv.lock
```

All scripts are run via `uv run python src/<script>.py` (see Makefile).

## Makefile Commands

All commands operate at the level of a `<City, Date Range>` pair:

```bash
# Set once and make all targets use them:
export CITY=chicago
export START_DATE=2024-01-01
export END_DATE=2024-03-31

make fetch        # BigQuery → /data/taveesh/scribe/raw/
make interpolate  # per-period interpolated grids
make aggregate    # per-period hex-aggregated GeoParquet
make cluster      # per-period cluster GeoJSON
make evaluate     # pairwise ARI stability score

make all          # interpolate + aggregate + cluster + evaluate
```

Key overridable variables (with defaults):

| Variable | Default | Options |
|---|---|---|
| `GRANULARITY` | `week` | `day`, `week`, `month` |
| `METHOD` | `idw` | `idw`, `loess`, `kde` |
| `RESOLUTION` | `8` | any H3 resolution |
| `N_CLUSTERS` | `auto` | integer or `auto` |
| `DISTANCE` | `Euclidean` | Manhattan, Euclidean, Chebyshev, Cosine, Correlation, Canberra |

For a development run with BigQuery sampling: `SAMPLE_PCT=1 make fetch ...`

BigQuery auth: set `GOOGLE_APPLICATION_CREDENTIALS=/path/to/key.json`.

## Architecture

```
src/
  cities.py            city name → bounding box / boundary polygon / UTM EPSG
  data_utils.py        load_parquet / save_parquet / to_geodataframe helpers
  geo_utils.py         H3 hexagon tessellation (h3 v4 API), grid generation
  interpolation.py     Interpolator base + IDW / LOESS_Interpolator / STBKR_Interpolator
  clustering.py        SpatialClustering base + SKATERClustering
  fetch_data.py        BigQuery ingestion script (CLI)
  run_interpolation.py interpolation script for a single period (CLI)
  aggregate.py         hex overlay + aggregation for a single period (CLI)
  run_clustering.py    SKATER clustering for a single period (CLI)
  evaluate.py          pairwise ARI stability across all periods (CLI)
  run_all.py           orchestrator that iterates over sub-periods and calls stage scripts
  seed_cities.py       one-time script to populate data/cities.geojson from OSM Nominatim

data/
  cities.geojson       committed; city boundary polygons keyed by city_key
```

### Key design decisions

- `fetch_data.py` queries the **full date range** once and caches it. Sub-period slicing is done locally in `run_interpolation.py` to avoid redundant BigQuery queries.
- `run_all.py` is the orchestrator called by Makefile; it splits the date range into sub-periods by granularity and invokes the single-period scripts for each.
- `aggregate.py` reuses `SpatialClustering.aggregate()` for consistent percentile computation across both interpolated-grid and clustering paths.
- `SKATERClustering` is called with `use_hexagons=False` in `run_clustering.py` because tessellation has already been done in `aggregate.py`.
- H3 v4 API is used throughout (`import h3`, not `from h3 import h3`). `cell_to_boundary()` returns `(lat, lng)` — all boundary-to-Polygon conversions flip to `(lng, lat)` for Shapely.

## Adding a New City

```bash
# Edit src/seed_cities.py and add the city to the CITIES list, then:
uv run python src/seed_cities.py
```

This fetches the city boundary polygon from OSM Nominatim and appends it to `data/cities.geojson`. Commit the updated GeoJSON.

## Storage Layout

```
/data/taveesh/scribe/
  raw/           {city}_{start}_{end}.parquet
  interpolated/  {city}_{period_start}_{period_end}_{method}.parquet
  aggregated/    {city}_{period_start}_{period_end}_{method}_res{N}.parquet
  output/        {city}_{period_start}_{period_end}_{method}_res{N}_clusters.geojson
                 {city}_{start}_{end}_{granularity}_{method}_stability.json
```
