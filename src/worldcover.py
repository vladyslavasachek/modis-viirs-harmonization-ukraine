"""Attach ESA WorldCover 2021 (v200) classes to detections and class shares to 0.1-deg cells.

Uses the COG internal overview at factor 8 (~67 m pixels) read over HTTP from the public AWS bucket,
cached as .npy in data/raw/worldcover/. Detection class = overview pixel containing the point.
Cell shares = fraction of overview pixels per class inside each 0.1-deg cell (Ukraine bbox).
Approximation: overview pixels are resampled from the 10 m map, so small patches (<~70 m) are smoothed out.
"""
import os

import numpy as np
import pandas as pd
import rasterio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "raw", "worldcover")
INTERIM = os.path.join(ROOT, "data", "interim")
PROC = os.path.join(ROOT, "data", "processed")
os.makedirs(CACHE, exist_ok=True)
URL = "/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"
F = 8
N = 36000 // F  # overview pixels per 3-deg tile
RES = 3.0 / N

WC_CLASSES = {10: "tree", 20: "shrub", 30: "grass", 40: "crop", 50: "built", 60: "bare",
              70: "snow", 80: "water", 90: "wetland", 95: "mangrove", 100: "moss"}


def tile_name(lat0, lon0):
    return f"N{lat0:02d}E{lon0:03d}"


def get_tile(lat0, lon0):
    path = os.path.join(CACHE, f"{tile_name(lat0, lon0)}_ov{F}.npy")
    if os.path.exists(path):
        return np.load(path)
    try:
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
            with rasterio.open(URL.format(tile=tile_name(lat0, lon0))) as s:
                a = s.read(1, out_shape=(N, N))
    except Exception as e:  # tile may not exist (sea)
        print("missing tile", tile_name(lat0, lon0), e)
        a = np.zeros((N, N), dtype=np.uint8)
    np.save(path, a)
    return a


def main():
    det = pd.read_parquet(os.path.join(INTERIM, "detections_all.parquet"))
    lat0 = (np.floor(det["lat"] / 3) * 3).astype(int)
    lon0 = (np.floor(det["lon"] / 3) * 3).astype(int)
    det["wc_code"] = np.uint8(0)
    tiles = sorted(set(zip(lat0, lon0)) | {(a, b) for a in range(42, 54, 3) for b in range(21, 42, 3)})
    cell_rows = []
    for (la, lo) in tiles:
        a = get_tile(la, lo)
        m = ((lat0 == la) & (lon0 == lo)).values
        if m.any():
            r = np.clip(((la + 3 - det.loc[m, "lat"].values) / RES).astype(int), 0, N - 1)
            c = np.clip(((det.loc[m, "lon"].values - lo) / RES).astype(int), 0, N - 1)
            det.loc[m, "wc_code"] = a[r, c]
        # cell shares: 0.1 deg = 150 overview px (3/0.1 = 30 cells per tile side)
        b = a.reshape(30, N // 30, 30, N // 30).transpose(0, 2, 1, 3).reshape(30, 30, -1)
        for code, name in WC_CLASSES.items():
            share = (b == code).mean(axis=2)
            for i in range(30):
                for j in range(30):
                    if share[i, j] > 0:
                        cell_lat = int(round((la + 3 - 0.1 * (i + 1)) * 10))
                        cell_lon = int(round((lo + 0.1 * j) * 10))
                        cell_rows.append((cell_lat * 1000 + cell_lon, name, share[i, j]))
        print("tile", tile_name(la, lo), "points", int(m.sum()))
    det["landcover"] = det["wc_code"].map(WC_CLASSES).fillna("nodata")
    det.to_parquet(os.path.join(INTERIM, "detections_all.parquet"), index=False)
    cs = pd.DataFrame(cell_rows, columns=["cell01", "lc", "share"])
    cs = cs.groupby(["cell01", "lc"])["share"].sum().unstack(fill_value=0).reset_index()
    cs.to_parquet(os.path.join(PROC, "cell01_landcover_shares.parquet"), index=False)
    print(det.groupby("sensor_product")["landcover"].value_counts(normalize=True).unstack().round(3).to_string())


if __name__ == "__main__":
    main()
