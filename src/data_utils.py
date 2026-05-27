import os
import pandas as pd
import geopandas as gpd


def load_parquet(path: str) -> pd.DataFrame:
    return pd.read_parquet(path)


def save_parquet(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_parquet(path, index=False)


def to_geodataframe(
    df: pd.DataFrame,
    lon_col: str,
    lat_col: str,
    crs: str = "EPSG:4326",
) -> gpd.GeoDataFrame:
    geometry = gpd.points_from_xy(df[lon_col], df[lat_col])
    return gpd.GeoDataFrame(df, geometry=geometry, crs=crs)
