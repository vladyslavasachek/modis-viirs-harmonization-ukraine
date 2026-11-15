"""Build compact app tables at 0.1 deg x week (B1 ratios are cell-independent, so they apply at any grid).

Outputs (data/processed/app/):
- weekly_cells.parquet : cell01, year, woy, month, source, h_all, h_excl, snpp_*, n20_*, modis_all, modis_excl
- cells.parquet        : cell01, lat, lon (centroid), oblast, in_cez, in_kyiv_airshed
- boundaries.geojson   : simplified oblasts + CEZ
"""
import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

sys.path.insert(0, os.path.dirname(__file__))
from harmonize import PROC, ROOT, load_det
from etl import CONFLICT_ZONES

OUT = os.path.join(PROC, "app")
os.makedirs(OUT, exist_ok=True)
OPEN_ENDED = [o for o, (a, b) in CONFLICT_ZONES.items() if b is None]
KYIV_AIRSHED = ["Kyiv", "Kyiv Oblast", "Zhytomyr Oblast", "Chernihiv Oblast"]


def main():
    d = load_det()
    d["excl"] = d["conflict_zone"] | d["oblast"].isin(OPEN_ENDED)
    P = np.load(os.path.join(PROC, "harmonization_params.npz"))
    R_ta, R_t, n20_ratio = P["R_ta"][0].reshape(24), P["R_t"][0].reshape(24), float(P["n20_ratio"])
    # source per (year, woy) from the Stage 2 table (uniform within a week)
    hc = pd.read_parquet(os.path.join(PROC, "harmonized_cells.parquet"), columns=["year", "woy", "source", "month"])
    src = hc.groupby(["year", "woy"]).agg(source=("source", "max"), month=("month", "min")).reset_index()

    keys = ["cell01", "year", "woy", "dn"]
    parts = {}
    ta = d["satellite"].isin(["Terra", "Aqua"])
    for nm, m in (("m_all", ta), ("m_excl", ta & ~d["excl"]), ("mT_all", d["satellite"] == "Terra"),
                  ("mT_excl", (d["satellite"] == "Terra") & ~d["excl"]),
                  ("snpp_all", d["satellite"] == "SNPP"), ("snpp_excl", (d["satellite"] == "SNPP") & ~d["excl"]),
                  ("n20_all", d["satellite"] == "NOAA20"), ("n20_excl", (d["satellite"] == "NOAA20") & ~d["excl"])):
        parts[nm] = d[m].groupby(keys).size()
    t = pd.DataFrame(parts).fillna(0).reset_index().merge(src, on=["year", "woy"], how="inner")
    strat = (t["month"].values - 1) * 2 + t["dn"].values
    for v in ("all", "excl"):
        t[f"h_{v}"] = np.select(
            [t["source"] == 0, t["source"] == 1, t["source"] == 2, t["source"] == 3],
            [t[f"snpp_{v}"], t[f"n20_{v}"] * n20_ratio, t[f"m_{v}"] * R_ta[strat], t[f"mT_{v}"] * R_t[strat]])
    num = ["m_all", "m_excl", "snpp_all", "snpp_excl", "n20_all", "n20_excl", "h_all", "h_excl"]
    w = t.groupby(["cell01", "year", "woy"]).agg(**{c: (c, "sum") for c in num},
                                                 month=("month", "first"), source=("source", "first")).reset_index()
    w = w[w[["h_all", "m_all", "snpp_all", "n20_all"]].sum(axis=1) > 0]
    for c in num[:6]:
        w[c] = w[c].astype("int16")
    for c in ("h_all", "h_excl"):
        w[c] = w[c].astype("float32")
    w = w.rename(columns={"m_all": "modis_all", "m_excl": "modis_excl"})
    w[["year", "woy", "month", "source"]] = w[["year", "woy", "month", "source"]].astype("int16")
    w.to_parquet(os.path.join(OUT, "weekly_cells.parquet"), index=False, compression="zstd")

    # cells + AOI membership (by cell centroid)
    cells = pd.DataFrame({"cell01": np.unique(w["cell01"])})
    cells["lat"] = (cells["cell01"] // 1000) / 10 + 0.05
    cells["lon"] = (cells["cell01"] % 1000) / 10 + 0.05
    adm1 = gpd.read_file(os.path.join(ROOT, "data", "raw", "boundaries", "UKR_ADM1.geojson"))[["shapeName", "geometry"]]
    pts = gpd.GeoDataFrame(cells, geometry=gpd.points_from_xy(cells["lon"], cells["lat"]), crs=4326)
    j = gpd.sjoin(pts, adm1, how="left", predicate="within")
    j = j[~j.index.duplicated()]
    # cells whose centroid falls just outside (coast/border): nearest oblast by detections
    mode_obl = (d.groupby(["cell01", "oblast"]).size().reset_index(name="k").sort_values("k")
                .drop_duplicates("cell01", keep="last").set_index("cell01")["oblast"])
    cells["oblast"] = j["shapeName"].values
    cells["oblast"] = cells["oblast"].fillna(cells["cell01"].map(mode_obl))
    cez = gpd.read_file(os.path.join(ROOT, "data", "raw", "cez", "cez.geojson")).geometry.iloc[0]
    cells["in_cez"] = shapely.contains_xy(cez, cells["lon"].values, cells["lat"].values)
    cells["in_kyiv_airshed"] = cells["oblast"].isin(KYIV_AIRSHED)
    cells.to_parquet(os.path.join(OUT, "cells.parquet"), index=False)

    b = adm1.copy()
    b["geometry"] = b.geometry.simplify(0.01)
    b["kind"] = "oblast"
    b = pd.concat([b, gpd.GeoDataFrame({"shapeName": ["Chornobyl Exclusion Zone"], "kind": ["cez"]},
                                       geometry=[cez.simplify(0.003)], crs=4326)])
    b.to_file(os.path.join(OUT, "boundaries.geojson"), driver="GeoJSON")
    print("weekly rows", len(w), "cells", len(cells), "CEZ cells", int(cells["in_cez"].sum()),
          {f: round(os.path.getsize(os.path.join(OUT, f)) / 1e6, 2) for f in os.listdir(OUT)})
    # ---------------- compact detection table with harmonization weights
    # weight = contribution of this detection to the harmonized (VIIRS-equivalent) value:
    # S-NPP weeks: S-NPP 1; NOAA-20-calibrated weeks: NOAA-20 ratio; MODIS-model weeks: B1 ratio(month, dn)
    # (Terra-only weeks: Terra ratio); all other detections 0. Sum of weights over any area/period = harmonized value.
    dd = d.merge(src, on=["year", "woy"], how="inner")
    mi = (dd["month"].values - 1) * 2 + dd["dn"].values
    sat = dd["satellite"].values
    s = dd["source"].values
    wgt = np.select([(s == 0) & (sat == "SNPP"), (s == 1) & (sat == "NOAA20"),
                     (s == 2) & np.isin(sat, ["Terra", "Aqua"]), (s == 3) & (sat == "Terra")],
                    [1.0, n20_ratio, R_ta[mi], R_t[mi]], 0.0)
    conf = dd["confidence"].astype(str)
    det = pd.DataFrame({
        "sensor": dd["satellite"].map({"Terra": "MODIS Terra", "Aqua": "MODIS Aqua", "SNPP": "VIIRS S-NPP",
                                       "NOAA20": "VIIRS NOAA-20", "NOAA21": "VIIRS NOAA-21"}).astype("category"),
        "datetime_utc": dd["datetime_utc"].dt.tz_localize(None).astype("datetime64[s]"),
        "lat": dd["lat"].astype("float32"), "lon": dd["lon"].astype("float32"),
        "frp": dd["frp"].round(1).astype("float32"),
        "confidence": conf.astype("category"),
        "daynight": dd["daynight"].astype("category"),
        "processing": dd["processing"].astype("category"),
        "conflict_flag": dd["excl"].values,
        "oblast": dd["oblast"].astype("category"),
        "landcover": dd["landcover"].astype("category"),
        "cell01": dd["cell01"].astype("int32"),
        "weight": wgt.astype("float32"),
    }).sort_values(["cell01", "datetime_utc"])
    det.to_parquet(os.path.join(OUT, "detections.parquet"), index=False, compression="zstd",
                   compression_level=9, row_group_size=200_000)
    print("detections", len(det), "MB", round(os.path.getsize(os.path.join(OUT, "detections.parquet")) / 1e6, 1),
          "weight-sum vs weekly h_all:", round(float(det["weight"].sum()), 1), round(float(w["h_all"].sum()), 1))
    # consistency check vs the 0.5-deg Stage 2 table
    hc2 = pd.read_parquet(os.path.join(PROC, "harmonized_cells.parquet"))
    a = w.groupby("year")[["h_all", "h_excl"]].sum()
    b2 = hc2.groupby("year")[["h_all", "h_excl"]].sum()
    print("max rel diff vs 0.5deg (annual):", float(((a - b2).abs() / b2).max().max()))


if __name__ == "__main__":
    main()
