import numpy as np
import pandas as pd
import geopandas as gpd
from joblib import Parallel, delayed
from tqdm import tqdm
from sklearn.neighbors import KernelDensity, NearestNeighbors
from geo_utils import assign_hex_ids, generate_grid_points, generate_grid_over_data, containment_filter
from loess.loess_2d import loess_2d


class Interpolator:
    """Base class for spatial interpolation methods.

    Parameters
    ----------
    data : pd.DataFrame
        Input measurement data.
    params : dict
        Keys: coord_columns (list[str]), h3_resolution (int), agg_fns (dict),
        grid_spacing (float), metric (str), epsg (int).
    shapefile : str, optional
        Path to city boundary shapefile. When provided, interpolation is done
        over a region-wide grid derived from the shapefile boundary.
    test_data : pd.DataFrame, optional
        Held-out data for validation. When provided, predictions are made at
        test-data locations rather than a synthetic grid.
    interpolate_to_grid : bool
        If True and shapefile is set, produce a region-wide grid output.
    """

    def __init__(
        self,
        data: pd.DataFrame,
        params: dict,
        shapefile: str = None,
        test_data: pd.DataFrame = None,
        interpolate_to_grid: bool = True,
    ):
        self.params = params
        self.coord_columns = params["coord_columns"]
        self.h3_resolution = params["h3_resolution"]
        self.agg_fns = params["agg_fns"]
        self.grid_spacing = params["grid_spacing"]
        self.metric = params["metric"]
        self.epsg = params["epsg"]

        lng, lat = data[self.coord_columns[0]], data[self.coord_columns[1]]
        self.data = gpd.GeoDataFrame(data, geometry=gpd.points_from_xy(lng, lat), crs="EPSG:4326")
        self.tessellation = None
        self.shapefile = shapefile
        self.interpolate_to_grid = interpolate_to_grid

        self.test_data = None
        if test_data is not None:
            tlng = test_data[self.coord_columns[0]]
            tlat = test_data[self.coord_columns[1]]
            self.test_data = gpd.GeoDataFrame(
                test_data,
                geometry=gpd.points_from_xy(tlng, tlat),
                crs="EPSG:4326",
            ).to_crs(self.epsg)

    def fit(self, n_jobs: int = -1) -> pd.DataFrame:
        """Fit the interpolation model.

        Returns a DataFrame with columns grid_x, grid_y, z_ground, z_smooth, w_smooth.
        When interpolate_to_grid is True or a shapefile is set, returns a single
        region-wide grid. Otherwise, returns per-hex-cell results.
        """
        self.data, self.tessellation = assign_hex_ids(self.data, self.h3_resolution)

        if self.shapefile is not None or self.interpolate_to_grid:
            grid_x, grid_y, z_ground, z_smooth, w_smooth = self.smooth_regionwide()
            return pd.DataFrame(
                {
                    "grid_x": grid_x,
                    "grid_y": grid_y,
                    "z_ground": z_ground,
                    "z_smooth": z_smooth,
                    "w_smooth": w_smooth,
                }
            )

        def process_cell(cell):
            cell_data = self.data[self.data["hex_id"] == cell]
            grid_x, grid_y, grid_z, z_smooth, w_smooth, coeff = self.smooth(cell)
            sample_size = len(cell_data)
            if z_smooth is None:
                agg_results = {k: np.nan for k in self.agg_fns}
                z_smooth = []
                w_smooth = []
            else:
                agg_results = {k: v(z_smooth) for k, v in self.agg_fns.items()}
                z_smooth = list(z_smooth)
                w_smooth = list(w_smooth)
            return {
                "h3_cell": cell,
                "sample_size": sample_size,
                "z_smooth": z_smooth,
                "grid_x": grid_x,
                "grid_y": grid_y,
                "grid_z": grid_z,
                "coeff": coeff,
                "w_smooth": w_smooth,
                **agg_results,
            }

        cells = list(self.data["hex_id"].unique())
        results = Parallel(n_jobs=n_jobs)(
            delayed(process_cell)(cell) for cell in tqdm(cells, desc="Smoothing", leave=False)
        )
        return pd.DataFrame(results)


class STBKR_Interpolator(Interpolator):
    """Adaptive KDE smoothing with k-NN bandwidth selection.

    Parameters
    ----------
    c : float
        Bandwidth scaling constant.
    k : int
        Number of nearest neighbours used to compute local bandwidth.
    """

    def __init__(
        self,
        data: pd.DataFrame,
        params: dict,
        c: float,
        k: int,
        shapefile: str = None,
        test_data: pd.DataFrame = None,
        interpolate_to_grid: bool = True,
    ):
        super().__init__(data, params, shapefile, test_data, interpolate_to_grid)
        self.c = c
        self.k = k

    def model_knn_distances(self, x: np.ndarray, y: np.ndarray, k: int) -> NearestNeighbors:
        k = min(k, len(x))
        knn = NearestNeighbors(n_neighbors=k)
        knn.fit(np.vstack((x, y)).T)
        return knn

    def kde_wrapper(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        xnew: np.ndarray,
        ynew: np.ndarray,
        c: float,
        k: int,
    ) -> tuple[np.ndarray | None, np.ndarray | None]:
        if len(x) == 0:
            return None, None
        knn = self.model_knn_distances(x, y, k)
        new_coords = np.vstack([xnew, ynew]).T
        distances, indices = knn.kneighbors(new_coords, return_distance=True)

        mean_distances = np.mean(distances, axis=1)
        bandwidths = 1e-12 + c * (mean_distances ** 2)

        z_smooth = np.zeros(len(xnew))
        for i in range(len(xnew)):
            bw = bandwidths[i]
            kde = KernelDensity(bandwidth=bw, kernel="gaussian")
            x_nb = x[indices[i]]
            y_nb = y[indices[i]]
            coords = np.vstack((x_nb, y_nb)).T
            kde.fit(coords)
            w = np.exp(kde.score_samples(coords))
            z_smooth[i] = np.sum(w * z[indices[i]]) / np.sum(w)

        w_smooth = np.ones(len(z_smooth))
        return z_smooth, w_smooth

    def smooth(self, h3_cell: str):
        cell_data = self.data[self.data["hex_id"] == h3_cell].to_crs(self.epsg)
        x = cell_data.geometry.x.values
        y = cell_data.geometry.y.values
        z = cell_data[self.metric].values

        if self.test_data is not None:
            td_within = containment_filter(h3_cell, self.test_data)
            td_within = td_within.to_crs(self.epsg)
            grid_x = td_within.geometry.x.values
            grid_y = td_within.geometry.y.values
            grid_z = td_within[self.metric].values
        else:
            grid_x, grid_y = generate_grid_points(h3_cell, self.epsg, self.grid_spacing)
            grid_z = np.full(len(grid_x), np.nan)

        z_smooth, w_smooth = self.kde_wrapper(x, y, z, xnew=grid_x, ynew=grid_y, c=self.c, k=self.k)
        if z_smooth is None:
            return grid_x, grid_y, grid_z, None, None, None

        coeff = np.full(z_smooth.shape, np.nan)
        return grid_x, grid_y, grid_z, z_smooth, w_smooth, coeff

    def smooth_regionwide(self):
        data = self.data.copy().to_crs(self.epsg)
        x = data.geometry.x.values
        y = data.geometry.y.values
        z = data[self.metric].values

        if not self.interpolate_to_grid and self.test_data is not None:
            xpt = self.test_data.geometry.x.values
            ypt = self.test_data.geometry.y.values
            zpt = self.test_data[self.metric].values
        else:
            xpt, ypt = generate_grid_over_data(self.epsg, self.grid_spacing, self.shapefile)
            zpt = np.full(len(xpt), np.nan)

        z_smooth, w_smooth = self.kde_wrapper(x, y, z, xnew=xpt, ynew=ypt, c=self.c, k=self.k)
        return xpt, ypt, zpt, z_smooth, w_smooth


class LOESS_Interpolator(Interpolator):
    """2D LOESS smoothing.

    Parameters
    ----------
    span : float
        Fraction of data points used in each local regression.
    degree : int
        Polynomial degree (1 or 2).
    """

    def __init__(
        self,
        data: pd.DataFrame,
        params: dict,
        span: float,
        degree: int,
        shapefile: str = None,
        test_data: pd.DataFrame = None,
        interpolate_to_grid: bool = True,
    ):
        super().__init__(data, params, shapefile, test_data, interpolate_to_grid)
        self.span = span
        self.degree = degree

    def smooth(self, h3_cell: str):
        cell_data = self.data[self.data["hex_id"] == h3_cell].to_crs(self.epsg)
        x = cell_data.geometry.x.values
        y = cell_data.geometry.y.values
        z = cell_data[self.metric].values

        if self.test_data is not None:
            td_within = containment_filter(h3_cell, self.test_data)
            td_within = td_within.to_crs(self.epsg)
            grid_x = td_within.geometry.x.values
            grid_y = td_within.geometry.y.values
            grid_z = td_within[self.metric].values
        else:
            grid_x, grid_y = generate_grid_points(h3_cell, self.epsg, self.grid_spacing)
            grid_z = np.full(len(grid_x), np.nan)

        try:
            z_smooth, w_smooth = loess_2d(
                x, y, z, xnew=grid_x, ynew=grid_y, frac=self.span, degree=self.degree
            )
        except np.linalg.LinAlgError:
            return grid_x, grid_y, grid_z, None, None, None

        return grid_x, grid_y, grid_z, z_smooth, w_smooth, None

    def smooth_regionwide(self):
        data = self.data.copy().to_crs(self.epsg)
        x = data.geometry.x.values
        y = data.geometry.y.values
        z = data[self.metric].values

        if not self.interpolate_to_grid and self.test_data is not None:
            xpt = self.test_data.geometry.x.values
            ypt = self.test_data.geometry.y.values
            zpt = self.test_data[self.metric].values
        else:
            xpt, ypt = generate_grid_over_data(self.epsg, self.grid_spacing, self.shapefile)
            zpt = np.full(len(xpt), np.nan)

        try:
            z_smooth, w_smooth = loess_2d(
                x, y, z, xnew=xpt, ynew=ypt, frac=self.span, degree=self.degree
            )
        except np.linalg.LinAlgError:
            return xpt, ypt, zpt, None, None

        return xpt, ypt, zpt, z_smooth, w_smooth


class IDW(Interpolator):
    """Inverse Distance Weighting interpolation.

    Parameters
    ----------
    p : float
        Power parameter. Higher values give more weight to nearby points.
    """

    def __init__(
        self,
        data: pd.DataFrame,
        params: dict,
        p: float = 2.0,
        shapefile: str = None,
        test_data: pd.DataFrame = None,
        interpolate_to_grid: bool = True,
    ):
        super().__init__(data, params, shapefile, test_data, interpolate_to_grid)
        self.p = p

    def predict(
        self,
        x0: np.ndarray,
        y0: np.ndarray,
        z0: np.ndarray,
        x1: np.ndarray,
        y1: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        predictions = []
        for xi, yi in zip(x1, y1):
            dist = np.sqrt((x0 - xi) ** 2 + (y0 - yi) ** 2)
            w = 1.0 / (dist + 1e-12) ** self.p
            w /= w.sum()
            predictions.append(np.dot(w, z0))
        return np.array(predictions), np.full(len(x1), np.nan)

    def smooth(self, h3_cell: str):
        cell_data = self.data[self.data["hex_id"] == h3_cell].to_crs(self.epsg)
        x = cell_data.geometry.x.values
        y = cell_data.geometry.y.values
        z = cell_data[self.metric].values

        if self.test_data is not None:
            td_within = containment_filter(h3_cell, self.test_data)
            td_within = td_within.to_crs(self.epsg)
            grid_x = td_within.geometry.x.values
            grid_y = td_within.geometry.y.values
            grid_z = td_within[self.metric].values
        else:
            grid_x, grid_y = generate_grid_points(h3_cell, self.epsg, self.grid_spacing)
            grid_z = np.full(len(grid_x), np.nan)

        z_smooth, w_smooth = self.predict(x, y, z, grid_x, grid_y)
        coeff = np.full(z_smooth.shape, np.nan)
        return grid_x, grid_y, grid_z, z_smooth, w_smooth, coeff

    def smooth_regionwide(self):
        data = self.data.copy().to_crs(self.epsg)
        x = data.geometry.x.values
        y = data.geometry.y.values
        z = data[self.metric].values

        if not self.interpolate_to_grid and self.test_data is not None:
            xpt = self.test_data.geometry.x.values
            ypt = self.test_data.geometry.y.values
            zpt = self.test_data[self.metric].values
        else:
            xpt, ypt = generate_grid_over_data(self.epsg, self.grid_spacing, self.shapefile)
            zpt = np.full(len(xpt), np.nan)

        z_smooth, w_smooth = self.predict(x, y, z, xpt, ypt)
        return xpt, ypt, zpt, z_smooth, w_smooth
