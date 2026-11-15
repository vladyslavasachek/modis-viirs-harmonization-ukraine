"""Stage 1 ETL: FIRMS CSVs -> clean detection table + base aggregate.

Outputs
- data/interim/detections_all.parquet   all filtered detections (incl. outside Ukraine; gitignored, large)
- data/processed/detections_ukr_viirs_snpp.parquet etc. are NOT written here (size); see aggregates
- data/processed/agg_cell_day.parquet   0.1 deg cell x UTC day x day/night x sensor (Ukraine only)
- data/processed/qc_*.csv               counts per filter step, sensitivity, version check
"""
from __future__ import annotations

import glob
import os

import geopandas as gpd
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
INTERIM = os.path.join(ROOT, "data", "interim")
PROC = os.path.join(ROOT, "data", "processed")
os.makedirs(INTERIM, exist_ok=True)
os.makedirs(PROC, exist_ok=True)

SENSOR_MAP = {  # FIRMS product code in file name -> (sensor, satellite-normalised)
    "M-C61": "MODIS",
    "SV-C2": "VIIRS_SNPP",
    "J1V-C2": "VIIRS_NOAA20",
    "J2V-C2": "VIIRS_NOAA21",
}

# Static-source mask grid (deg). ~0.004 deg ~ 375-450 m at 48N.
MASK_RES = 0.004
MASK_MIN_HITS = 3  # cell is "persistent static source" if >=3 type-2 detections (any sensor, standard data)

# Conflict-affected zones: oblast-level, date windows (approximation; documented in METHODS).
CONFLICT_START = pd.Timestamp("2022-02-24", tz="UTC")
CONFLICT_ZONES = {
    # northern front: fighting / occupation Feb 24 - early April 2022 (withdrawal completed ~Apr 3-7)
    "Kyiv Oblast": ("2022-02-24", "2022-04-07"),
    "Kyiv": ("2022-02-24", "2022-04-07"),
    "Chernihiv Oblast": ("2022-02-24", "2022-04-07"),
    "Sumy Oblast": ("2022-02-24", "2022-04-07"),
    "Zhytomyr Oblast": ("2022-02-24", "2022-04-07"),
    # east and south: open-ended from invasion start
    "Kharkiv Oblast": ("2022-02-24", None),
    "Donetsk Oblast": ("2022-02-24", None),
    "Luhansk Oblast": ("2022-02-24", None),
    "Zaporizhia Oblast": ("2022-02-24", None),
    "Kherson Oblast": ("2022-02-24", None),
    "Mykolaiv Oblast": ("2022-02-24", None),
    "Autonomous Republic of Crimea": ("2022-02-24", None),
    "Sevastopol": ("2022-02-24", None),
}


def load_raw() -> pd.DataFrame:
    frames = []
    for path in sorted(glob.glob(os.path.join(RAW, "firms", "*", "fire_*.csv"))):
        name = os.path.basename(path)
        processing = "standard" if name.startswith("fire_archive") else "nrt"
        code = next(k for k in SENSOR_MAP if f"_{k}_" in name)
        df = pd.read_csv(path, dtype={"confidence": str, "version": str, "acq_time": int})
        df["sensor_product"] = SENSOR_MAP[code]
        df["processing"] = processing
        df["source_file"] = name
        if "type" not in df:
            df["type"] = np.nan
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)

    df["sensor"] = np.where(df["instrument"] == "MODIS", "MODIS", "VIIRS")
    df["satellite"] = df["satellite"].replace({"N20": "NOAA20", "N21": "NOAA21", "SNPP": "SNPP"})
    t = df["acq_time"].astype(int)
    df["datetime_utc"] = pd.to_datetime(df["acq_date"]) + pd.to_timedelta(t // 100, "h") + pd.to_timedelta(t % 100, "min")
    df["datetime_utc"] = df["datetime_utc"].dt.tz_localize("UTC")
    df["date"] = df["datetime_utc"].dt.floor("D").dt.tz_localize(None)
    df = df.rename(columns={"latitude": "lat", "longitude": "lon"})
    df["type"] = df["type"].astype("float").astype("Int8")
    return df[
        ["sensor", "satellite", "sensor_product", "datetime_utc", "date", "daynight", "lat", "lon", "frp",
         "scan", "track", "brightness", "confidence", "type", "processing", "version", "source_file"]
    ]


def static_mask_keys(df: pd.DataFrame) -> set:
    s = df[(df["processing"] == "standard") & (df["type"] == 2)]
    k = list(zip(np.floor(s["lat"] / MASK_RES).astype(int), np.floor(s["lon"] / MASK_RES).astype(int)))
    vc = pd.Series(k).value_counts()
    return set(vc[vc >= MASK_MIN_HITS].index)


def add_admin(df: pd.DataFrame) -> pd.DataFrame:
    adm1 = gpd.read_file(os.path.join(RAW, "boundaries", "UKR_ADM1.geojson"))[["shapeName", "geometry"]]
    pts = gpd.GeoDataFrame(df[["lat"]].copy(), geometry=gpd.points_from_xy(df["lon"], df["lat"]), crs=4326)
    j = gpd.sjoin(pts, adm1, how="left", predicate="within")
    j = j[~j.index.duplicated()]
    df["oblast"] = j["shapeName"].values
    df["in_ukraine"] = df["oblast"].notna()
    return df


def add_conflict_flag(df: pd.DataFrame) -> pd.DataFrame:
    flag = np.zeros(len(df), dtype=bool)
    for obl, (a, b) in CONFLICT_ZONES.items():
        m = (df["oblast"] == obl) & (df["date"] >= pd.Timestamp(a))
        if b is not None:
            m &= df["date"] <= pd.Timestamp(b)
        flag |= m.values
    df["conflict_zone"] = flag
    return df


def confidence_ok(df: pd.DataFrame, modis_min=30, viirs_set=("n", "h")) -> pd.Series:
    is_m = df["sensor"] == "MODIS"
    mconf = pd.to_numeric(df["confidence"].where(is_m), errors="coerce")
    return (is_m & (mconf >= modis_min)) | (~is_m & df["confidence"].isin(viirs_set))


def main():
    df = load_raw()
    qc = []

    def log(step, d):
        g = d.groupby(["sensor_product", "processing"]).agg(rows=("lat", "size"), first=("date", "min"), last=("date", "max"))
        g["step"] = step
        qc.append(g.reset_index())

    log("0_raw", df)

    # Static sources: standard -> by type; NRT -> by mask built from standard type 2
    keys = static_mask_keys(df)
    k = list(zip(np.floor(df["lat"] / MASK_RES).astype(int), np.floor(df["lon"] / MASK_RES).astype(int)))
    df["static_mask"] = pd.Series(k, index=df.index).isin(keys)
    is_std = df["processing"] == "standard"
    veg = (is_std & (df["type"] == 0)) | (~is_std & ~df["static_mask"])
    print(f"static-mask cells: {len(keys)}; NRT rows masked: {int((~is_std & df['static_mask']).sum())}")
    # how well does the mask reproduce type on standard data? (QC)
    std = df[is_std]
    tp = ((std["type"] == 2) & std["static_mask"]).sum()
    print(f"mask on standard: recall of type2 {tp / (std['type'] == 2).sum():.3f}; "
          f"type0 wrongly masked {((std['type'] == 0) & std['static_mask']).mean() / (std['type'] == 0).mean():.4f}")
    df = df[veg].copy()
    log("1_vegetation", df)

    # Sensitivity of confidence thresholds (on vegetation detections)
    sens = []
    for mmin in (0, 30, 50, 80):
        for vset, vlab in ((("l", "n", "h"), "lnh"), (("n", "h"), "nh"), (("h",), "h")):
            ok = confidence_ok(df, mmin, vset)
            for sp, n in df[ok].groupby("sensor_product").size().items():
                sens.append({"modis_min": mmin, "viirs_set": vlab, "sensor_product": sp, "rows": n})
    pd.DataFrame(sens).to_csv(os.path.join(PROC, "qc_confidence_sensitivity.csv"), index=False)

    df = df[confidence_ok(df)].copy()
    log("2_confidence", df)

    df = add_admin(df)
    df = add_conflict_flag(df)
    log("3_in_ukraine", df[df["in_ukraine"]])

    df["cell01"] = (np.floor(df["lat"] * 10).astype(int) * 1000 + np.floor(df["lon"] * 10).astype(int))
    df.to_parquet(os.path.join(INTERIM, "detections_all.parquet"), index=False)

    pd.concat(qc).to_csv(os.path.join(PROC, "qc_filter_steps.csv"), index=False)
    print(pd.concat(qc).to_string())
    print("conflict-flagged (Ukraine):", int(df.loc[df.in_ukraine, "conflict_zone"].sum()))


if __name__ == "__main__":
    main()
