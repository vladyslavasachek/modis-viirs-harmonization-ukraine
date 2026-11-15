"""Data access + AOI computations for the Streamlit app (uses only data/processed/)."""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import shapely
import streamlit as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from harmonize import QS, SOURCES, aoi_components, aoi_scale, summarize_components  # noqa: E402

PROC = os.path.join(ROOT, "data", "processed")
APP = os.path.join(PROC, "app")

DOMAIN_LABELS = {
    "excl": "Ukraine without conflict-affected oblasts, all years",
    "all": "All of Ukraine, incl. conflict-affected oblasts",
}
CONFLICT_NOTE = (
    "**Ukraine without conflict-affected oblasts, all years** removes Kharkiv, Donetsk, Luhansk, Zaporizhzhia, "
    "Kherson, Mykolaiv oblasts, Crimea and Sevastopol in *every* year (so years stay comparable), and Kyiv, "
    "Chernihiv, Sumy and Zhytomyr oblasts only during 24 Feb – 7 Apr 2022. These are approximate oblast-level "
    "zones of hostilities since Feb 2022; we do **not** claim that any individual detection is caused by "
    "military activity.")
PRESETS_FIXED = ["Ukraine", "Kyiv airshed (Kyiv city + Kyiv, Zhytomyr, Chernihiv oblasts)",
                 "Chornobyl Exclusion Zone", "Custom polygon (draw on map)"]


@st.cache_data
def load_weekly() -> pd.DataFrame:
    return pd.read_parquet(os.path.join(APP, "weekly_cells.parquet"))


@st.cache_data
def load_cells() -> pd.DataFrame:
    return pd.read_parquet(os.path.join(APP, "cells.parquet"))


@st.cache_data
def load_boundaries() -> dict:
    with open(os.path.join(APP, "boundaries.geojson")) as f:
        return json.load(f)


@st.cache_data
def load_calib():
    ta = pd.read_csv(os.path.join(PROC, "interval_calibration_B1_TA.csv"))
    t = pd.read_csv(os.path.join(PROC, "interval_calibration_B1_T.csv"))
    n20q = pd.read_csv(os.path.join(PROC, "noaa20_ratio_rel_quantiles.csv"), index_col=0).iloc[:, 0].to_dict()
    n20_ratio = float(np.load(os.path.join(PROC, "harmonization_params.npz"))["n20_ratio"])
    return ta, t, n20q, n20_ratio


@st.cache_data
def all_weeks() -> pd.DataFrame:
    w = load_weekly().groupby(["year", "woy"]).agg(month=("month", "min"), source=("source", "max")).reset_index()
    w["week_start"] = pd.to_datetime(w["year"].astype(str) + "-01-01") + pd.to_timedelta(w["woy"] * 7, "D")
    return w


def oblast_names() -> list[str]:
    return sorted(load_cells()["oblast"].dropna().unique())


def aoi_cells(preset: str, polygon=None) -> np.ndarray:
    c = load_cells()
    if preset == "Ukraine":
        m = np.ones(len(c), bool)
    elif preset.startswith("Kyiv airshed"):
        m = c["in_kyiv_airshed"].values
    elif preset == "Chornobyl Exclusion Zone":
        m = c["in_cez"].values
    elif preset.startswith("Custom"):
        if polygon is None:
            return np.array([], dtype=int)
        m = shapely.contains_xy(polygon, c["lon"].values, c["lat"].values)
    else:
        m = (c["oblast"] == preset).values
    return c.loc[m, "cell01"].values


@st.cache_data(show_spinner=False)
def aoi_weekly(cells_key: tuple, domain: str) -> pd.DataFrame:
    """Weekly harmonized series for the AOI with 80/95% intervals and source labels."""
    cells = np.array(cells_key)
    w = load_weekly()
    t = w[w["cell01"].isin(cells)]
    calib_ta, calib_t, n20q, n20_ratio = load_calib()
    weeks = all_weeks()
    if len(t):
        comp = aoi_components(t, domain, n20_ratio)
    else:
        comp = pd.DataFrame(columns=["year", "woy", "obs", "n20", "mod_ta", "mod_t", "point"])
    comp = weeks[["year", "woy", "month", "source"]].merge(
        comp.drop(columns=[c for c in ("month", "source", "has_ta", "has_t") if c in comp]),
        on=["year", "woy"], how="left").fillna(0)
    comp["has_ta"] = comp["source"] == 2
    comp["has_t"] = comp["source"] == 3
    scale = aoi_scale(int(np.ceil(len(cells) / 25)))
    out = summarize_components(comp, ["year", "woy"], calib_ta, calib_t, n20q, scale)
    out = out.merge(weeks[["year", "woy", "month", "source", "week_start"]], on=["year", "woy"])
    out["source_label"] = out["source"].map(SOURCES)
    raw = t.groupby(["year", "woy"])[[f"modis_{domain}", f"snpp_{domain}", f"n20_{domain}"]].sum()
    raw.columns = ["raw_modis", "raw_snpp", "raw_noaa20"]
    out = out.merge(raw.reset_index(), on=["year", "woy"], how="left").fillna(
        {"raw_modis": 0, "raw_snpp": 0, "raw_noaa20": 0})
    out.attrs["scale"] = scale
    return out


def aggregate(weekly: pd.DataFrame, keys) -> pd.DataFrame:
    """Aggregate a weekly AOI series to months/years with recomputed intervals."""
    calib_ta, calib_t, n20q, _ = load_calib()
    scale = weekly.attrs.get("scale", "ukr")
    a = summarize_components(weekly, keys, calib_ta, calib_t, n20q, scale)
    raw = weekly.groupby(keys)[["raw_modis", "raw_snpp", "raw_noaa20"]].sum().reset_index()
    return a.merge(raw, on=keys)


def fmt_interval(df: pd.DataFrame) -> pd.Series:
    def f(r):
        if r["share_estimated"] > 0.001 if pd.notna(r["share_estimated"]) else False:
            return f"{r['point']:,.0f} (80%: {r['lo80']:,.0f}–{r['hi80']:,.0f}; 95%: {r['lo95']:,.0f}–{r['hi95']:,.0f})"
        return f"{r['point']:,.0f} (observed)"
    return df.apply(f, axis=1)


# ----------------------------------------------------------------------------- detections / KPIs / trends
UNIT = "VIIRS-equivalent hot spots"
PNG_CONFIG = {"toImageButtonOptions": {"format": "png", "scale": 2}, "displaylogo": False}


@st.cache_data
def load_detections() -> pd.DataFrame:
    """Decode the compact detection table (see src/build_app_data.py)."""
    d = pd.read_parquet(os.path.join(APP, "detections.parquet"))
    d["lat"] = (d["lat"] / 1e4).astype("float32")
    d["lon"] = (d["lon"] / 1e4).astype("float32")
    d["frp"] = (d["frp"] / 10).astype("float32")
    d["datetime_utc"] = pd.to_datetime(d["datetime_utc"].astype("int64") * 60, unit="s")
    d["cell01"] = (np.floor(d["lat"].astype(float) * 10).astype("int32") * 1000
                   + np.floor(d["lon"].astype(float) * 10).astype("int32"))
    return d


def haversine_km(lat0, lon0, lat, lon):
    lat0, lon0, lat, lon = map(np.radians, (lat0, lon0, lat, lon))
    a = np.sin((lat - lat0) / 2) ** 2 + np.cos(lat0) * np.cos(lat) * np.sin((lon - lon0) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def circle_cells(lat0, lon0, radius_km) -> np.ndarray:
    c = load_cells()
    return c.loc[haversine_km(lat0, lon0, c["lat"].values, c["lon"].values) <= radius_km, "cell01"].values


def aoi_detections(aoi: dict, domain: str) -> pd.DataFrame:
    det = load_detections()
    if aoi["kind"] == "circle":
        box = det[(det["lat"].between(aoi["lat"] - aoi["r"] / 100, aoi["lat"] + aoi["r"] / 100))]
        dist = haversine_km(aoi["lat"], aoi["lon"], box["lat"].values.astype(float), box["lon"].values.astype(float))
        out = box[dist <= aoi["r"]]
    elif aoi["kind"] == "polygon":
        b = shapely.bounds(aoi["polygon"])
        box = det[det["lat"].between(b[1], b[3]) & det["lon"].between(b[0], b[2])]
        out = box[shapely.contains_xy(aoi["polygon"], box["lon"].values, box["lat"].values)]
    else:
        out = det[det["cell01"].isin(aoi["cells"])]
    if domain == "excl":
        out = out[~out["conflict_flag"]]
    return out


def month_to_date_kpi(det: pd.DataFrame, n20q: dict, base_years=range(2001, 2026)) -> dict | None:
    """Harmonized activity for the current (latest, partial) month-to-date vs the same calendar window in
    earlier years. Value = sum of harmonization weights."""
    allmax = load_detections()["datetime_utc"].max()
    end = pd.Timestamp(allmax.date())
    start = end.replace(day=1)
    d = det.assign(y=det["datetime_utc"].dt.year, m=det["datetime_utc"].dt.month, dd=det["datetime_utc"].dt.day)
    win = d[(d["m"] == end.month) & (d["dd"] <= end.day)]
    per_year = win.groupby("y")["weight"].sum().reindex(list(base_years) + [end.year], fill_value=0.0)
    cur = float(per_year.loc[end.year])
    cur_n20 = float(win[(win["y"] == end.year) & (win["sensor"] == "VIIRS NOAA-20")]["weight"].sum())
    hist = per_year.loc[list(base_years)].values
    med = float(np.median(hist))
    rng = np.random.default_rng(0)
    boots = [np.median(rng.choice(hist, len(hist))) for _ in range(1000)]
    pct = float((np.sum(hist < cur) + 0.5 * np.sum(hist == cur)) / len(hist) * 100)
    return dict(start=start, end=end, value=cur, value_lo=cur - cur_n20 * (1 - n20q["lo95"]),
                value_hi=cur + cur_n20 * (n20q["hi95"] - 1), median=med,
                median_lo=float(np.quantile(boots, .025)), median_hi=float(np.quantile(boots, .975)),
                percentile=pct, per_year=per_year)


def sen_mk(years: np.ndarray, values: np.ndarray) -> dict:
    from scipy.stats import kendalltau, theilslopes
    s = theilslopes(values, years, alpha=0.95)
    tau, p = kendalltau(years, values)
    return dict(slope=s.slope, lo=s.low_slope, hi=s.high_slope, tau=tau, p=p,
                intercept=s.intercept, rel=s.slope / np.median(values) * 100 if np.median(values) > 0 else np.nan)


OPEN_ENDED = ["Kharkiv Oblast", "Donetsk Oblast", "Luhansk Oblast", "Zaporizhia Oblast", "Kherson Oblast",
              "Mykolaiv Oblast", "Autonomous Republic of Crimea", "Sevastopol"]


def annual_from_weights(det: pd.DataFrame) -> pd.DataFrame:
    """Annual harmonized value (sum of weights) and raw single-sensor counts from a detection subset."""
    y = det["datetime_utc"].dt.year
    out = pd.DataFrame({"harmonized": det["weight"].groupby(y).sum(),
                        "aqua": (det["sensor"] == "MODIS Aqua").groupby(y).sum(),
                        "snpp": (det["sensor"] == "VIIRS S-NPP").groupby(y).sum()})
    out.index.name = "year"
    return out


def constant_area(det_all: pd.DataFrame) -> pd.DataFrame:
    """Constant-area trend domain: conflict-affected oblasts with open-ended windows removed in ALL years;
    the north Feb-Apr 2022 window is kept (excluding it would shrink the area in 2022 only)."""
    return det_all[~det_all["oblast"].isin(OPEN_ENDED)]


def trend_robustness(det_all: pd.DataFrame) -> pd.DataFrame:
    """Sen's slope + Mann-Kendall on harmonized vs raw single-sensor series, complete years, fixed periods.
    det_all: AOI detections WITHOUT domain filtering. Raw series are in their own sensor's units: compare sign
    and % of median per year, not levels."""
    rows = []
    groups = (("Constant area: without conflict-affected oblasts (all years)", annual_from_weights(constant_area(det_all)),
               [("Harmonized", "harmonized", 2001, 2025), ("Raw Aqua-MODIS only", "aqua", 2003, 2025),
                ("Raw VIIRS S-NPP only", "snpp", 2012, 2025)]),
              ("Constant area, pre-invasion only", annual_from_weights(constant_area(det_all)),
               [("Harmonized", "harmonized", 2001, 2021)]),
              ("All of Ukraine / AOI, pre-invasion", annual_from_weights(det_all),
               [("Harmonized", "harmonized", 2001, 2021), ("Raw Aqua-MODIS only", "aqua", 2003, 2021),
                ("Raw VIIRS S-NPP only", "snpp", 2012, 2021)]),
              ("All of Ukraine / AOI, incl. conflict-affected oblasts", annual_from_weights(det_all),
               [("Harmonized", "harmonized", 2001, 2025)]))
    for dom, ann, specs in groups:
        for label, col, y0, y1 in specs:
            yrs = np.arange(y0, y1 + 1)
            v = ann[col].reindex(yrs, fill_value=0).values.astype(float)
            if (v > 0).sum() < 6:
                continue
            t = sen_mk(yrs, v)
            rows.append({"domain": dom, "series": label, "period": f"{y0}–{y1}", "Sen's slope / yr": t["slope"],
                         "95% CI": f"{t['lo']:+,.0f} to {t['hi']:+,.0f}", "% of median / yr": t["rel"],
                         "Mann–Kendall τ": t["tau"], "p": t["p"]})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- burn types / events
TYPE_ORDER = ["crop residue / field", "grass / abandoned land", "forest", "peat / wetland smouldering",
              "other / uncertain"]
TYPE_COLORS = {"crop residue / field": "#E69F00", "grass / abandoned land": "#56B4E9", "forest": "#009E73",
               "peat / wetland smouldering": "#D55E00", "other / uncertain": "#999999"}   # Okabe-Ito (CVD-safe)


@st.cache_data
def load_events() -> pd.DataFrame:
    e = pd.read_parquet(os.path.join(APP, "events.parquet"))
    e["cell01"] = (np.floor(e["lat"].astype(float) * 10).astype(int) * 1000
                   + np.floor(e["lon"].astype(float) * 10).astype(int))
    return e


def aoi_events(aoi: dict, domain: str, d0, d1) -> pd.DataFrame:
    e = load_events()
    e = e[(e["start"] >= d0) & (e["start"] < d1 + pd.Timedelta(days=1))]
    if aoi["kind"] == "circle":
        e = e[haversine_km(aoi["lat"], aoi["lon"], e["lat"].values.astype(float), e["lon"].values.astype(float)) <= aoi["r"]]
    elif aoi["kind"] == "polygon":
        e = e[shapely.contains_xy(aoi["polygon"], e["lon"].values, e["lat"].values)]
    else:
        e = e[e["cell01"].isin(aoi["cells"])]
    if domain == "excl":
        e = e[~e["oblast"].isin(OPEN_ENDED) & ~e["conflict_event"]]
    return e


def event_percentiles(ev: pd.Series) -> pd.DataFrame:
    """Percentile of an event's size/intensity/duration among all events (>=2 detections) in the same oblast."""
    e = load_events()
    ref = e[e["oblast"] == ev["oblast"]]
    rows = []
    for col, label, unit in (("n_det", "Detections (size)", ""), ("extent_km", "Spatial extent", " km"),
                             ("frp_sum", "Total FRP (energy proxy)", " MW"), ("frp_max", "Peak FRP", " MW"),
                             ("duration_days", "Duration", " days")):
        v = float(ev[col])
        r = ref[col].values.astype(float)
        pct = (np.sum(r < v) + 0.5 * np.sum(r == v)) / len(r) * 100
        rows.append({"metric": label, "this event": f"{v:,.1f}{unit}", "oblast median": f"{np.median(r):,.1f}{unit}",
                     "percentile in oblast": pct})
    return pd.DataFrame(rows), len(ref)


WV_LAYERS = ("Coastlines_15m,Reference_Features_15m(hidden),Reference_Labels_15m(hidden),"
             "VIIRS_SNPP_Thermal_Anomalies_375m_All,VIIRS_NOAA20_Thermal_Anomalies_375m_All,"
             "MODIS_Combined_Thermal_Anomalies_All,MODIS_Aqua_CorrectedReflectance_TrueColor,"
             "MODIS_Terra_CorrectedReflectance_TrueColor")


def worldview(lat, lon, date, pad=0.25):
    """NASA Worldview link (same layers as src/events_validation.py)."""
    return (f"https://worldview.earthdata.nasa.gov/?v={lon - pad:.3f},{lat - pad * 0.65:.3f},{lon + pad:.3f},"
            f"{lat + pad * 0.65:.3f}&l={WV_LAYERS}&t={date:%Y-%m-%d}-T12:00:00Z")


def weekly_from_detections(det: pd.DataFrame, value: str = "weight", sensor: str | None = None) -> pd.DataFrame:
    """Weekly series (year, woy, point) from detections: harmonized (sum of weights) or raw counts of one sensor."""
    d = det if sensor is None else det[det["sensor"] == sensor]
    t = d["datetime_utc"]
    g = pd.DataFrame({"year": t.dt.year, "woy": np.minimum((t.dt.dayofyear - 1) // 7, 51),
                      "v": d[value].astype(float) if sensor is None else 1.0})
    return g.groupby(["year", "woy"])["v"].sum().rename("point").reset_index()
