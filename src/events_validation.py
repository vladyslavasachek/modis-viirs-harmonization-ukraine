"""Stage 5 validation outputs + app tables.

1. Known events: CEZ April 2020; Kyiv-region (Vyshhorod) peat fires Aug-Oct 2015.
2. 40 events for visual checking in NASA Worldview: stratified random (8 per rule type) among events with >= 3
   detections, 2012-2025, not in conflict-affected zones. Column `my_label` left empty.
3. App tables: events.parquet reduced to events with >= 2 detections; per-detection burn type added to
   detections.parquet (S-NPP: event type; other sensors / pre-2012 MODIS: land-cover-only type, lower confidence).
"""
import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

sys.path.insert(0, os.path.dirname(__file__))
from events import LC_TO_TYPE, TYPES
from harmonize import INTERIM, PROC, ROOT

APP = os.path.join(PROC, "app")
# MODIS true colour as base (S-NPP imagery layers are flagged after the end of S-NPP delivery, Nov 2026);
# fire overlays: VIIRS S-NPP 375 m, VIIRS NOAA-20 375 m (from 2018), MODIS Terra/Aqua.
LAYERS = ("Coastlines_15m,Reference_Features_15m(hidden),Reference_Labels_15m(hidden),"
          "VIIRS_SNPP_Thermal_Anomalies_375m_All,VIIRS_NOAA20_Thermal_Anomalies_375m_All,"
          "MODIS_Combined_Thermal_Anomalies_All,MODIS_Aqua_CorrectedReflectance_TrueColor,"
          "MODIS_Terra_CorrectedReflectance_TrueColor")


def worldview(lat, lon, date, pad=0.25):
    return (f"https://worldview.earthdata.nasa.gov/?v={lon - pad:.3f},{lat - pad * 0.65:.3f},{lon + pad:.3f},"
            f"{lat + pad * 0.65:.3f}&l={LAYERS}&t={date:%Y-%m-%d}-T12:00:00Z")


def main():
    e = pd.read_parquet(os.path.join(INTERIM, "events_all.parquet"))

    # ---- 1. known events
    cez = gpd.read_file(os.path.join(ROOT, "data", "raw", "cez", "cez.geojson")).geometry.iloc[0]
    rows = []
    for name, m in (
        ("Chornobyl Exclusion Zone fires, Apr 2020",
         shapely.contains_xy(cez, e["lon"].values, e["lat"].values) & (e["start"] >= "2020-03-25") & (e["start"] < "2020-05-01")),
        ("Kyiv-region peat fires, Vyshhorod district, 15 Aug-31 Oct 2015",
         e["lat"].between(50.5, 51.0) & e["lon"].between(30.1, 30.9) & (e["start"] >= "2015-08-15") & (e["start"] < "2015-11-01")),
    ):
        s = e[m]
        for t, g in s.groupby("burn_type", observed=True):
            rows.append(dict(case=name, burn_type=t, events=len(g), detections=int(g["n_det"].sum()),
                             share_of_detections=g["n_det"].sum() / s["n_det"].sum()))
    known = pd.DataFrame(rows)
    known.to_csv(os.path.join(PROC, "events_known_cases.csv"), index=False)
    print(known.round(3).to_string(index=False))

    # ---- 2. 40-event Worldview sample
    pool = e[(e["n_det"] >= 3) & (e["year"] <= 2025) & (~e["conflict_event"])]
    rng = np.random.default_rng(2026)
    picks = []
    for t in TYPES:
        g = pool[pool["burn_type"] == t]
        picks.append(g.iloc[rng.choice(len(g), min(8, len(g)), replace=False)])
    s = pd.concat(picks).sort_values("start")
    out = pd.DataFrame({
        "event_id": s["event_id"].values, "start_utc": s["start"].dt.strftime("%Y-%m-%d %H:%M").values,
        "end_utc": s["end"].dt.strftime("%Y-%m-%d %H:%M").values, "lat": s["lat"].round(4).values,
        "lon": s["lon"].round(4).values, "oblast": s["oblast"].astype(str).values, "n_detections": s["n_det"].values,
        "duration_days": s["duration_days"].round(2).values, "frp_mean_MW": s["frp_mean"].round(1).values,
        "night_fraction": s["night_frac"].round(2).values,
        "share_crop": s["lc_crop"].round(2).values, "share_grass": s["lc_grass"].round(2).values,
        "share_tree": s["lc_tree"].round(2).values, "share_wetland": s["lc_wetland"].round(2).values,
        "rule_label": s["burn_type"].astype(str).values, "label_confidence": s["type_confidence"].astype(str).values,
        "gmm_majority_label": s["gmm_label"].astype(str).values,
        "worldview_during": [worldview(a, b, d) for a, b, d in zip(s["lat"], s["lon"], s["start"])],
        "worldview_after_burn_scar": [worldview(a, b, d + pd.Timedelta(days=3)) for a, b, d in
                                      zip(s["lat"], s["lon"], s["end"])],
        "my_label": "", "my_notes": ""})
    vpath = os.path.join(ROOT, "validation_events_worldview.csv")
    if os.path.exists(vpath):   # never overwrite the manually labelled file
        print("validation CSV exists - not overwritten (event ids are deterministic)")
    else:
        out.to_csv(vpath, index=False)
        print("validation sample:", out["rule_label"].value_counts().to_dict())

    # ---- 3. app tables
    small = e[e["n_det"] >= 2].drop(columns=["gmm_cluster"])
    small.to_parquet(os.path.join(APP, "events.parquet"), index=False, compression="zstd", compression_level=19)
    print("events (>=2 detections)", len(small), "MB", round(os.path.getsize(os.path.join(APP, "events.parquet")) / 1e6, 2))
    single = e[e["n_det"] == 1]
    pd.Series({"events_total": len(e), "events_ge2": len(small), "singletons": len(single)}).to_csv(
        os.path.join(PROC, "events_counts.csv"))

    # per-detection burn type
    mem = pd.read_parquet(os.path.join(INTERIM, "viirs_event_members.parquet"))
    mem = mem.merge(e[["event_id", "burn_type"]], on="event_id", how="left")
    mem["datetime_utc"] = mem["datetime_utc"].dt.tz_localize(None).astype("datetime64[s]")
    # same float32 rounding path as src/build_app_data.py
    mem["lat_i"] = np.round(mem["lat"].astype("float32").astype(float) * 1e4).astype("int32")
    mem["lon_i"] = np.round(mem["lon"].astype("float32").astype(float) * 1e4).astype("int32")
    mem["t_i"] = (mem["datetime_utc"].astype("int64") // 60).astype("int32")
    mem = mem.drop_duplicates(["t_i", "lat_i", "lon_i"])
    det = pd.read_parquet(os.path.join(APP, "detections.parquet"))
    if "burn_type" in det:
        det = det.drop(columns=[c for c in ("burn_type", "type_from_event", "event_id") if c in det])
    is_v = det["sensor"].isin(["VIIRS S-NPP", "VIIRS NOAA-20"]).values
    m = det.loc[is_v, ["datetime_utc", "lat", "lon"]].reset_index().rename(
        columns={"datetime_utc": "t_i", "lat": "lat_i", "lon": "lon_i"})
    m = m.merge(mem[["t_i", "lat_i", "lon_i", "burn_type", "event_id"]], on=["t_i", "lat_i", "lon_i"], how="left")
    m = m.drop_duplicates("index")
    bt = det["landcover"].astype(str).map(LC_TO_TYPE).fillna(TYPES[4]).values.astype(object)
    from_event = np.zeros(len(det), bool)
    ev_id = np.full(len(det), -1, dtype=np.int32)
    hit = m["burn_type"].notna().values
    bt[m.loc[hit, "index"].values] = m.loc[hit, "burn_type"].astype(str).values
    from_event[m.loc[hit, "index"].values] = True
    ev_id[m.loc[hit, "index"].values] = m.loc[hit, "event_id"].values
    det["burn_type"] = pd.Categorical(bt, categories=TYPES)
    det["type_from_event"] = from_event
    det.to_parquet(os.path.join(APP, "detections.parquet"), index=False, compression="zstd", compression_level=19,
                   row_group_size=1_000_000)
    snpp = det["sensor"] == "VIIRS S-NPP"
    print("S-NPP detections matched to events:", round(float(det.loc[snpp, "type_from_event"].mean()), 4),
          "detections MB", round(os.path.getsize(os.path.join(APP, "detections.parquet")) / 1e6, 1))
    w = det.groupby("burn_type", observed=True)["weight"].sum()
    print("harmonized weight by type:", (w / w.sum()).round(3).to_dict())


if __name__ == "__main__":
    main()
