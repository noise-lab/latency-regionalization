##
## Latency Regionalization Pipeline
##
## Usage:
##   make fetch        CITY=chicago START_DATE=2024-01-01 END_DATE=2024-03-31
##   make interpolate  CITY=chicago START_DATE=2024-01-01 END_DATE=2024-03-31
##   make aggregate    CITY=chicago START_DATE=2024-01-01 END_DATE=2024-03-31
##   make cluster      CITY=chicago START_DATE=2024-01-01 END_DATE=2024-03-31
##   make evaluate     CITY=chicago START_DATE=2024-01-01 END_DATE=2024-03-31
##   make all          CITY=chicago START_DATE=2024-01-01 END_DATE=2024-03-31
##
## Optional overrides (shown with defaults):
##   GRANULARITY=week       # day | week | month
##   METHOD=idw             # idw | loess | kde
##   RESOLUTION=8           # H3 resolution
##   MIN_SAMPLES=5          # minimum grid points per hex cell
##   N_CLUSTERS=auto        # integer or auto
##   DISTANCE=Euclidean     # Manhattan | Euclidean | Chebyshev | Cosine | Correlation | Canberra
##   N_BOOTSTRAPS=100
##   GRID_SPACING=500       # metres in projected CRS
##   IDW_P=2                # IDW power parameter
##   LOESS_SPAN=0.3
##   LOESS_DEGREE=2
##   KDE_C=0.5
##   KDE_K=10
##

CITY         ?= chicago
START_DATE   ?= 2024-01-01
END_DATE     ?= 2024-12-31
GRANULARITY  ?= week
METHOD       ?= idw
RESOLUTION   ?= 8
MIN_SAMPLES  ?= 5
N_CLUSTERS   ?= auto
DISTANCE     ?= Euclidean
N_BOOTSTRAPS ?= 100
GRID_SPACING ?= 500
IDW_P        ?= 2
LOESS_SPAN   ?= 0.3
LOESS_DEGREE ?= 2
KDE_C        ?= 0.5
KDE_K        ?= 10

RUN = uv run python src

# Arguments forwarded to run_all.py
_COMMON = \
	--city $(CITY) \
	--start $(START_DATE) \
	--end $(END_DATE) \
	--granularity $(GRANULARITY) \
	--method $(METHOD) \
	--resolution $(RESOLUTION)

_INTERP_ARGS = \
	--grid-spacing $(GRID_SPACING) \
	--p $(IDW_P) \
	--span $(LOESS_SPAN) \
	--degree $(LOESS_DEGREE) \
	--c $(KDE_C) \
	--k $(KDE_K)

_AGG_ARGS = --min-samples $(MIN_SAMPLES)

_CLUSTER_ARGS = \
	--n-clusters $(N_CLUSTERS) \
	--distance-metric $(DISTANCE) \
	--n-bootstraps $(N_BOOTSTRAPS)

# ── Targets ────────────────────────────────────────────────────────────────────

.PHONY: fetch interpolate aggregate cluster evaluate all setup

## Fetch raw M-Lab data from BigQuery (full date range, cached as Parquet).
fetch:
	$(RUN)/fetch_data.py \
		--city $(CITY) \
		--start $(START_DATE) \
		--end $(END_DATE)

## Interpolate raw measurements to a regular grid for each sub-period.
interpolate:
	$(RUN)/run_all.py $(_COMMON) --stage interpolate $(_INTERP_ARGS)

## Overlay H3 hexagons and compute per-cell latency aggregates.
aggregate:
	$(RUN)/run_all.py $(_COMMON) --stage aggregate $(_AGG_ARGS)

## Run SKATER spatial clustering on hex-aggregated data.
cluster:
	$(RUN)/run_all.py $(_COMMON) --stage cluster $(_CLUSTER_ARGS)

## Evaluate clustering stability across sub-periods (pairwise ARI).
evaluate:
	$(RUN)/evaluate.py \
		--city $(CITY) \
		--start $(START_DATE) \
		--end $(END_DATE) \
		--granularity $(GRANULARITY) \
		--method $(METHOD) \
		--resolution $(RESOLUTION)

## Run interpolate → aggregate → cluster → evaluate in sequence.
all: interpolate aggregate cluster evaluate

## Install dependencies with uv.
setup:
	uv sync
