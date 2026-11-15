"""ETL for extra FIRMS regions used only for the Kyiv upwind fire load (e.g. data/raw/firms_belarus/).

Any folder data/raw/firms_<region>/ is picked up automatically (no code change needed). Same filters as Ukraine:
standard processing type == 0; NRT: static-source mask built from that region's standard type-2 detections
(0.004 deg, >= 3 hits); VIIRS confidence n/h; MODIS confidence >= 30. WorldCover class attached; burn types from
the same event clustering + rules as Ukraine (S-NPP stream, then NOAA-20 after S-NPP ends).
Output: data/interim/extra_detections.parquet
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from etl import MASK_MIN_HITS, MASK_RES, SENSOR_MAP, confidence_ok
from events import LC_TO_TYPE, TYPES, cluster, features, rule_types
from harmonize import INTERIM, ROOT
from worldcover import RES, WC_CLASSES, get_tile


def load_region(folder):
    frames = []
    for path in sorted(glob.glob(os.path.join(folder, "*", "fire_*.csv"))):
        name = os.path.basename(path)
        code = next(k for k in SENSOR_MAP if f"_{k}_" in name)
        df = pd.read_csv(path, dtype={"confidence": str, "version": str, "acq_time": int})
        df["processing"] = "standard" if name.startswith("fire_archive") else "nrt"
        df["sensor_product"] = SENSOR_MAP[code]
        if "type" not in df:
            df["type"] = np.nan
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df["sensor"] = np.where(df["instrument"] == "MODIS", "MODIS", "VIIRS")
    df["satellite"] = df["satellite"].replace({"N20": "NOAA20", "N21": "NOAA21"})
    t = df["acq_time"].astype(int)
    df["datetime_utc"] = (pd.to_datetime(df["acq_date"]) + pd.to_timedelta(t // 100, "h")
                          + pd.to_timedelta(t % 100, "min")).dt.tz_localize("UTC")
    df["date"] = df["datetime_utc"].dt.floor("D").dt.tz_localize(None)
    return df.rename(columns={"latitude": "lat", "longitude": "lon"})


def main():
    out = []
    for folder in sorted(glob.glob(os.path.join(ROOT, "data", "raw", "firms_*"))):
        region = os.path.basename(folder).replace("firms_", "").capitalize()
        df = load_region(folder)
        n0 = len(df)
        s = df[(df["processing"] == "standard") & (df["type"] == 2)]
        k = pd.Series(list(zip(np.floor(s["lat"] / MASK_RES).astype(int), np.floor(s["lon"] / MASK_RES).astype(int))))
        vc = k.value_counts()
        keys = set(vc[vc >= MASK_MIN_HITS].index)
        kk = pd.Series(list(zip(np.floor(df["lat"] / MASK_RES).astype(int), np.floor(df["lon"] / MASK_RES).astype(int))))
        static = kk.isin(keys).values
        std = df["processing"] == "standard"
        df = df[((std & (df["type"] == 0)) | (~std & ~static)).values]
        df = df[confidence_ok(df)].copy()
        # WorldCover
        lat0 = (np.floor(df["lat"] / 3) * 3).astype(int)
        lon0 = (np.floor(df["lon"] / 3) * 3).astype(int)
        df["wc_code"] = 0
        for la, lo in sorted(set(zip(lat0, lon0))):
            a = get_tile(la, lo)
            m = ((lat0 == la) & (lon0 == lo)).values
            r = np.clip(((la + 3 - df.loc[m, "lat"].values) / RES).astype(int), 0, a.shape[0] - 1)
            c = np.clip(((df.loc[m, "lon"].values - lo) / RES).astype(int), 0, a.shape[1] - 1)
            df.loc[m, "wc_code"] = a[r, c]
        df["landcover"] = df["wc_code"].map(WC_CLASSES).fillna("nodata")
        df["oblast"] = region
        df["conflict_zone"] = False
        # events + types on the S-NPP stream (NOAA-20 after S-NPP ends)
        snpp_end = df.loc[df["satellite"] == "SNPP", "date"].max()
        v = df[(df["satellite"] == "SNPP") | ((df["satellite"] == "NOAA20") & (df["date"] > snpp_end))]
        v = v.sort_values("datetime_utc").reset_index()
        v["event_id"] = cluster(v)
        v["recur"] = False   # no MODIS history downloaded for extra regions (feature unused by the rules)
        e = rule_types(features(v))
        v = v.merge(e[["event_id", "burn_type"]], on="event_id")
        df["burn_type"] = df["landcover"].map(LC_TO_TYPE).fillna(TYPES[4])
        df["type_from_event"] = False
        df.loc[v["index"].values, "burn_type"] = v["burn_type"].values
        df.loc[v["index"].values, "type_from_event"] = True
        df["region"] = region
        print(f"{region}: raw {n0:,} -> filtered {len(df):,}; events {e['event_id'].nunique():,}; "
              f"{df['date'].min():%Y-%m-%d}..{df['date'].max():%Y-%m-%d}", df["satellite"].value_counts().to_dict())
        print(df.groupby("burn_type")["frp"].agg(["size", "sum"]).round(0).to_string())
        out.append(df[["region", "satellite", "processing", "datetime_utc", "date", "daynight", "lat", "lon", "frp",
                       "landcover", "burn_type", "type_from_event"]])
    if out:
        pd.concat(out).to_parquet(os.path.join(INTERIM, "extra_detections.parquet"), index=False)


if __name__ == "__main__":
    main()
