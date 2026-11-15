"""PEATMAP (Xu et al. 2018, CATENA 160; dataset doi:10.5518/252, CC BY 4.0) peatland polygons for Ukraine/Belarus."""
import os

import geopandas as gpd
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHP = os.path.join(ROOT, "data", "raw", "peatmap", "Europe", "Other_European_Peatland.shp")


def on_peat(lat, lon) -> np.ndarray:
    """True where the point falls inside a PEATMAP polygon (False everywhere if PEATMAP is missing)."""
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    if not os.path.exists(SHP):
        return np.zeros(len(lat), bool)
    g = gpd.read_file(SHP).to_crs(4326)
    g = g.cx[lon.min() - 0.1:lon.max() + 0.1, lat.min() - 0.1:lat.max() + 0.1][["geometry"]]
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon, lat), crs=4326)
    j = gpd.sjoin(pts, g, how="left", predicate="within")
    hit = j["index_right"].notna().groupby(level=0).any()
    return hit.reindex(range(len(lat)), fill_value=False).values
