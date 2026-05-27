import spopt
from sklearn.metrics import pairwise as skm
import libpysal
from geo_utils import assign_hex_ids
import numpy as np
import geopandas as gpd
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd


class SpatialClustering:
    """Base class for spatial clustering."""

    def __init__(
        self,
        data: gpd.GeoDataFrame,
        resolution: int,
        aggregates: list,
        metric: str,
        use_hexagons: bool,
        n_bootstraps: int = 1000,
    ):
        self.data = data
        self.metric = metric
        self.resolution = resolution
        self.aggregates = aggregates
        self.use_hexagons = use_hexagons
        self.n_bootstraps = n_bootstraps

    def aggregate(
        self, sample: gpd.GeoDataFrame, over: str, tessellation: gpd.GeoDataFrame
    ) -> gpd.GeoDataFrame:
        grouped = sample.groupby(over).agg(
            min=(self.metric, "min"),
            max=(self.metric, "max"),
            mean=(self.metric, "mean"),
            std=(self.metric, lambda x: np.std(x) if len(x) > 1 else 0),
            p10=(self.metric, lambda x: np.percentile(x, 10)),
            p25=(self.metric, lambda x: np.percentile(x, 25)),
            p50=(self.metric, lambda x: np.percentile(x, 50)),
            p75=(self.metric, lambda x: np.percentile(x, 75)),
            p90=(self.metric, lambda x: np.percentile(x, 90)),
            p95=(self.metric, lambda x: np.percentile(x, 95)),
            p975=(self.metric, lambda x: np.percentile(x, 97.5)),
            p99=(self.metric, lambda x: np.percentile(x, 99)),
            latency_reduction=(self.metric, lambda x: np.percentile(x, 90) - np.percentile(x, 10)),
            norm_latency_reduction=(
                self.metric,
                lambda x: (np.percentile(x, 90) - np.percentile(x, 10)) / np.percentile(x, 10),
            ),
            inequality_ratio=(self.metric, lambda x: np.percentile(x, 90) / np.percentile(x, 10)),
        ).reset_index()
        geom = tessellation[[over, "geometry"]].drop_duplicates()
        geometries = dict(zip(geom[over], geom["geometry"]))
        return gpd.GeoDataFrame(grouped, geometry=[geometries[h] for h in grouped[over]])

    @classmethod
    def calculate_jaccard_index(cls, a: set, b: set) -> float:
        union = a.union(b)
        if not union:
            return 0.0
        return len(a.intersection(b)) / len(union)

    def aggregate_over_tessellation(self, over: str = "hex_id") -> None:
        self.tessellation = self.aggregate(self.data, over, self.tessellation)

    def tessellate(
        self, sample: gpd.GeoDataFrame
    ) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
        if self.use_hexagons:
            sample, tessellation = assign_hex_ids(sample, self.resolution)
        else:
            tessellation = sample.copy()
        if tessellation.crs is None:
            tessellation = tessellation.set_crs("EPSG:4326")
        else:
            tessellation = tessellation.to_crs("EPSG:4326")
        return sample, tessellation


class SKATERClustering(SpatialClustering):
    """Spatial clustering using the SKATER algorithm."""

    def __init__(
        self,
        data: gpd.GeoDataFrame,
        resolution: int,
        aggregates: list,
        metric: str,
        use_hexagons: bool,
        n_bootstraps: int = 100,
    ):
        super().__init__(data, resolution, aggregates, metric, use_hexagons, n_bootstraps)
        self.metric_func = {
            "Manhattan": skm.manhattan_distances,
            "Euclidean": skm.euclidean_distances,
            "Chebyshev": lambda x, y: skm.pairwise_distances(x, y, metric="chebyshev"),
            "Cosine": skm.cosine_distances,
            "Correlation": lambda x, y: skm.pairwise_distances(x, y, metric="correlation"),
            "Canberra": lambda x, y: skm.pairwise_distances(x, y, metric="canberra"),
        }

    def build_model(
        self,
        distance_metric: str,
        tessellation: gpd.GeoDataFrame,
        weights,
        n_clusters: int,
        floor: int,
    ) -> spopt.region.Skater:
        spanning_forest_kwds = dict(
            dissimilarity=self.metric_func[distance_metric],
            affinity=None,
            reduction=np.sum,
            center=np.mean,
        )
        return spopt.region.Skater(
            tessellation,
            weights,
            self.aggregates,
            n_clusters=n_clusters,
            floor=floor,
            trace=False,
            islands="ignore",
            spanning_forest_kwds=spanning_forest_kwds,
        )

    def cluster(
        self,
        n_clusters: int,
        floor: int,
        distance_metric: str,
        over: str = "hex_id",
    ) -> gpd.GeoDataFrame:
        data, tessellation = self.tessellate(self.data)
        tessellation = self.aggregate(data, over, tessellation)
        self.weights = libpysal.weights.Queen.from_dataframe(tessellation, use_index=False)
        model = self.build_model(distance_metric, tessellation, self.weights, n_clusters, floor)
        model.solve()
        tessellation["cluster"] = model.labels_
        tessellation["n_clusters"] = n_clusters
        tessellation["floor"] = floor
        tessellation["distance_metric"] = distance_metric
        tessellation["algorithm"] = "SKATER"
        return tessellation

    def get_optimal_clusters(
        self, distance_metric: str, over: str = "hex_id"
    ) -> tuple[int, int]:
        """Search K=1..10, stopping when consecutive cluster assignments are identical."""
        data, tessellation = self.tessellate(self.data)
        tessellation = self.aggregate(data, over, tessellation)
        num_units = tessellation[over].nunique()
        # floor = 10% of hex units; ensures clusters are large enough to be meaningful
        floor = num_units // 10
        prev = self.cluster(1, floor, distance_metric, over)["cluster"].values
        for K in range(2, 11):
            curr = self.cluster(K, floor, distance_metric, over)["cluster"].values
            if np.array_equal(prev, curr):
                return K - 1, floor
            prev = curr
        return 10, floor

    def plot(self, clus: gpd.GeoDataFrame) -> None:
        fig, ax = plt.subplots(1, 1)
        clus["cluster"] = clus["cluster"].astype(str)
        clus.plot(column="cluster", cmap="tab20", legend=False, edgecolor="black", linewidth=0.5, ax=ax)
        plt.axis("off")

    def calculate_bootstrap_cluster_iteration(
        self,
        iteration: int,
        n_clusters: int,
        floor: int,
        distance_metric: str,
        over: str = "hex_id",
    ) -> gpd.GeoDataFrame:
        sample = self.data.sample(frac=1, replace=True)
        sample, tessellation = self.tessellate(sample)
        tessellation = self.aggregate(sample, over, tessellation)
        weights = libpysal.weights.Queen.from_dataframe(tessellation, use_index=False)
        model = self.build_model(distance_metric, tessellation, weights, n_clusters, floor)
        model.solve()
        tessellation["cluster"] = model.labels_
        tessellation["n_clusters"] = n_clusters
        tessellation["floor"] = floor
        tessellation["distance_metric"] = distance_metric
        tessellation["algorithm"] = "SKATER"
        tessellation["bootstrap_id"] = iteration
        return tessellation

    def calculate_bootstrap_clusters(
        self,
        n_clusters: int,
        floor: int,
        distance_metric: str,
        over: str = "hex_id",
    ) -> pd.DataFrame:
        if self.data.crs is None:
            self.data = self.data.set_crs("EPSG:4326")
        else:
            self.data = self.data.to_crs("EPSG:4326")
        data = [
            self.calculate_bootstrap_cluster_iteration(i, n_clusters, floor, distance_metric, over)
            for i in tqdm(range(self.n_bootstraps), desc="Bootstrapping...")
        ]
        return pd.concat(data, axis=0)

    def calculate_jaccard_matrix(
        self, tessellation: pd.DataFrame, n_clusters: int, over: str = "hex_id"
    ) -> np.ndarray:
        jaccard_matrix = np.zeros((n_clusters, self.n_bootstraps, self.n_bootstraps))
        for cid in range(n_clusters):
            for i in range(self.n_bootstraps):
                for j in range(i, self.n_bootstraps):
                    a = set(tessellation[(tessellation["bootstrap_id"] == i) & (tessellation["cluster"] == cid)][over].values)
                    b = set(tessellation[(tessellation["bootstrap_id"] == j) & (tessellation["cluster"] == cid)][over].values)
                    jaccard_matrix[cid, i, j] = self.calculate_jaccard_index(a, b)
                    jaccard_matrix[cid, j, i] = jaccard_matrix[cid, i, j]
        return jaccard_matrix

    def plot_jaccard_matrix(self, jaccard_matrix: np.ndarray) -> None:
        plt.imshow(jaccard_matrix, cmap="viridis")
        plt.colorbar()
        plt.show()
