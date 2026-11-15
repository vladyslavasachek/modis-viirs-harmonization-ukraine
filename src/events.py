"""Stage 5: fire events by space-time clustering of VIIRS detections + rule-based burn types + GMM cross-check.

Clustering (ST-DBSCAN with min_samples=1, i.e. connected components of the space-time neighbour graph):
two detections are linked if they are <= EPS_KM apart in space AND <= EPS_DAYS apart in time.
- EPS_KM = 1.0: VIIRS I-band pixels are 375 m at nadir and up to ~800 m at the scan edge; 1 km links adjacent
  pixels (incl. diagonals and off-nadir growth) without merging separate fields.
- EPS_DAYS = 2: bridges one missed/clouded overpass (S-NPP sees Ukraine ~2x per day).
Input stream: VIIRS S-NPP (reference sensor, consistent over 2012-2026); after the end of S-NPP (2 Nov 2026)
NOAA-20 NRT. Clustering is done per calendar year (events crossing 31 Dec are split; negligible in winter).
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(__file__))
from harmonize import INTERIM, PROC
from peat import on_peat

EPS_KM, EPS_DAYS = 1.0, 2.0
OUT = os.path.join(PROC, "app")
TYPES = ["crop residue / field", "grass / abandoned land", "forest", "peat / wetland smouldering", "other / uncertain"]
LC_TO_TYPE = {"crop": TYPES[0], "grass": TYPES[1], "shrub": TYPES[1], "tree": TYPES[2], "wetland": TYPES[3]}


def load():
    d = pd.read_parquet(os.path.join(INTERIM, "detections_all.parquet"),
                        columns=["satellite", "datetime_utc", "date", "daynight", "lat", "lon", "frp", "in_ukraine",
                                 "conflict_zone", "oblast", "landcover", "processing"])
    d = d[d["in_ukraine"]]
    snpp_end = d.loc[d["satellite"] == "SNPP", "date"].max()
    v = d[(d["satellite"] == "SNPP") | ((d["satellite"] == "NOAA20") & (d["date"] > snpp_end))].copy()
    v = v.sort_values("datetime_utc").reset_index(drop=True)
    return d, v


def cluster(v: pd.DataFrame) -> np.ndarray:
    labels = np.full(len(v), -1, dtype=np.int64)
    nxt = 0
    t = (v["datetime_utc"] - pd.Timestamp("2012-01-01", tz="UTC")).dt.total_seconds().values / 86400
    x = v["lon"].values * 111.32 * np.cos(np.radians(v["lat"].values))
    y = v["lat"].values * 110.57
    yr = v["datetime_utc"].dt.year.values
    for Y in np.unique(yr):
        idx = np.where(yr == Y)[0]
        P = np.column_stack([x[idx] / EPS_KM, y[idx] / EPS_KM, t[idx] / EPS_DAYS])
        pairs = cKDTree(P).query_pairs(r=1.0, p=np.inf, output_type="ndarray")   # box: |dx|,|dy|<=1 km, |dt|<=2 d
        if len(pairs):
            ds = np.hypot(x[idx][pairs[:, 0]] - x[idx][pairs[:, 1]], y[idx][pairs[:, 0]] - y[idx][pairs[:, 1]])
            pairs = pairs[ds <= EPS_KM]                                            # circle in space
        n = len(idx)
        g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)) if len(pairs) else \
            coo_matrix((n, n))
        nc, lab = connected_components(g, directed=False)
        labels[idx] = lab + nxt
        nxt += nc
        print(f"{Y}: {n:,} detections -> {nc:,} events", flush=True)
    return labels


def recurrence(v: pd.DataFrame, d_all: pd.DataFrame) -> np.ndarray:
    """Per detection: did its ~1 km cell (0.01 deg) have any MODIS or VIIRS detection in the previous 5 years?
    MODIS history from 2001 makes this comparable for all events 2012+."""
    def key(lat, lon):
        return np.floor(lat * 100).astype(np.int64) * 100000 + np.floor(lon * 100).astype(np.int64)
    hist = pd.DataFrame({"k": key(d_all["lat"].values, d_all["lon"].values), "y": d_all["date"].dt.year.values})
    hist = hist.drop_duplicates()
    ky = pd.DataFrame({"k": key(v["lat"].values, v["lon"].values), "y": v["datetime_utc"].dt.year.values})
    out = np.zeros(len(v), bool)
    for lag in range(1, 6):
        h = hist.assign(y=hist["y"] + lag)
        out |= ky.merge(h.assign(hit=True), on=["k", "y"], how="left")["hit"].fillna(False).values.astype(bool)
    return out


def features(v: pd.DataFrame) -> pd.DataFrame:
    v = v.assign(night=(v["daynight"] == "N").astype(float), conflict=v["conflict_zone"].astype(float))
    if "on_peat" not in v:
        v["on_peat"] = on_peat(v["lat"].values, v["lon"].values)
    v["on_peat"] = v["on_peat"].astype(float)
    for lc in ("crop", "grass", "tree", "wetland", "built"):
        v[f"lc_{lc}"] = (v["landcover"] == lc).astype(float)
    v["lc_other"] = 1 - v[[f"lc_{c}" for c in ("crop", "grass", "tree", "wetland")]].sum(axis=1)
    x = v["lon"].values * 111.32 * np.cos(np.radians(v["lat"].values))
    y = v["lat"].values * 110.57
    v["x"], v["y"] = x, y
    g = v.groupby("event_id")
    e = g.agg(start=("datetime_utc", "min"), end=("datetime_utc", "max"), n_det=("frp", "size"),
              frp_mean=("frp", "mean"), frp_max=("frp", "max"), frp_sum=("frp", "sum"),
              night_frac=("night", "mean"), lat=("lat", "mean"), lon=("lon", "mean"),
              lc_crop=("lc_crop", "mean"), lc_grass=("lc_grass", "mean"), lc_tree=("lc_tree", "mean"),
              lc_wetland=("lc_wetland", "mean"), lc_built=("lc_built", "mean"), lc_other=("lc_other", "mean"),
              recur_frac=("recur", "mean"), conflict_frac=("conflict", "mean"), peat_frac=("on_peat", "mean"),
              xmin=("x", "min"), xmax=("x", "max"), ymin=("y", "min"), ymax=("y", "max"))
    e["oblast"] = g["oblast"].agg(lambda s: s.mode().iat[0] if s.notna().any() else None) if len(e) < 2000 else \
        v.groupby(["event_id", "oblast"]).size().reset_index(name="k").sort_values("k").drop_duplicates(
            "event_id", keep="last").set_index("event_id")["oblast"].reindex(e.index)
    e["duration_days"] = ((e["end"] - e["start"]).dt.total_seconds() / 86400).round(2)
    e["extent_km"] = np.hypot(e["xmax"] - e["xmin"], e["ymax"] - e["ymin"]).round(2)
    e = e.drop(columns=["xmin", "xmax", "ymin", "ymax"])
    e["doy"] = e["start"].dt.dayofyear
    e["year"] = e["start"].dt.year
    return e.reset_index()


def _rules(e, peat):
    lab = np.full(len(e), TYPES[4], dtype=object)
    forest = ~peat & (e["lc_tree"] >= 0.5)
    crop = ~peat & ~forest & (e["lc_crop"] >= 0.5) & (e["duration_days"] <= 3)
    grass = ~peat & ~forest & ~crop & (e["lc_grass"] >= 0.5)
    for m, t in ((peat, TYPES[3]), (forest, TYPES[2]), (crop, TYPES[0]), (grass, TYPES[1])):
        lab[np.asarray(m)] = t
    return lab


def rule_types(e: pd.DataFrame) -> pd.DataFrame:
    """Transparent rules (applied in this order):
    1 peat / wetland smouldering: (a) >= 50% of detections on mapped peat (PEATMAP, Xu et al. 2018), regardless of
      WorldCover - drained peatlands are often mapped as grass/cropland; OR (b) wetland share >= 0.3 AND
      duration >= 2 d AND night fraction >= 0.35 AND mean FRP <= 10 MW (long, night-persistent, low-intensity)
    2 forest: tree share >= 0.5
    3 crop residue / field: cropland share >= 0.5 AND duration <= 3 d
    4 grass / abandoned land: grass (+shrub) share >= 0.5
    5 other / uncertain: everything else (mixed land cover, long cropland events, built-up, ...)
    Confidence: high if the deciding land-cover share >= 0.8, medium if >= 0.6, else low (peat: always medium,
    because PEATMAP is a coarse compilation and WorldCover herbaceous wetland is only a proxy)."""
    lab = np.full(len(e), TYPES[4], dtype=object)
    share = np.zeros(len(e))
    peat_sig = (e["lc_wetland"] >= 0.3) & (e["duration_days"] >= 2) & (e["night_frac"] >= 0.35) & (e["frp_mean"] <= 10)
    peat_map = e["peat_frac"] >= 0.5 if "peat_frac" in e else False
    # label without the peat map (for before/after comparison)
    e["burn_type_wc_only"] = _rules(e, peat_sig)
    peat = peat_sig | peat_map
    forest = ~peat & (e["lc_tree"] >= 0.5)
    crop = ~peat & ~forest & (e["lc_crop"] >= 0.5) & (e["duration_days"] <= 3)
    grass = ~peat & ~forest & ~crop & (e["lc_grass"] >= 0.5)
    for m, t, col in ((peat, TYPES[3], "lc_wetland"), (forest, TYPES[2], "lc_tree"), (crop, TYPES[0], "lc_crop"),
                      (grass, TYPES[1], "lc_grass")):
        lab[m.values] = t
        share[m.values] = e.loc[m, col].values
    conf = np.where(share >= 0.8, "high", np.where(share >= 0.6, "medium", "low"))
    conf[np.asarray(peat)] = "medium"
    e["burn_type"] = lab
    e["type_confidence"] = conf
    return e


def gmm_crosscheck(e: pd.DataFrame, k: int = 5, seed: int = 0):
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    from sklearn.mixture import GaussianMixture
    from sklearn.preprocessing import StandardScaler
    cols = ["log_n", "log_dur", "log_frp", "night_frac", "lc_crop", "lc_grass", "lc_tree", "lc_wetland",
            "doy_sin", "doy_cos", "recur_frac"]
    X = pd.DataFrame({"log_n": np.log1p(e["n_det"]), "log_dur": np.log1p(e["duration_days"]),
                      "log_frp": np.log1p(e["frp_mean"]), "night_frac": e["night_frac"],
                      "lc_crop": e["lc_crop"], "lc_grass": e["lc_grass"], "lc_tree": e["lc_tree"],
                      "lc_wetland": e["lc_wetland"], "doy_sin": np.sin(2 * np.pi * e["doy"] / 365.25),
                      "doy_cos": np.cos(2 * np.pi * e["doy"] / 365.25), "recur_frac": e["recur_frac"]})[cols]
    Z = StandardScaler().fit_transform(X)
    rng = np.random.default_rng(seed)
    fit_idx = rng.choice(len(Z), min(len(Z), 150_000), replace=False)
    gm = GaussianMixture(k, covariance_type="full", random_state=seed, n_init=2).fit(Z[fit_idx])
    cl = gm.predict(Z)
    e["gmm_cluster"] = cl
    ct = pd.crosstab(e["gmm_cluster"], e["burn_type"])
    major = ct.idxmax(axis=1)
    e["gmm_label"] = e["gmm_cluster"].map(major)
    stats = dict(k=k, ARI=adjusted_rand_score(e["burn_type"], cl), NMI=normalized_mutual_info_score(e["burn_type"], cl),
                 agreement=float((e["gmm_label"] == e["burn_type"]).mean()))
    prof = X.assign(cluster=cl).groupby("cluster").mean()
    return e, ct, stats, prof


def main():
    t0 = time.time()
    d_all, v = load()
    v["event_id"] = cluster(v)
    v["recur"] = recurrence(v, d_all)
    print(f"clustered {len(v):,} detections in {time.time() - t0:.0f}s", flush=True)
    e = rule_types(features(v))
    e, ct, stats, prof = gmm_crosscheck(e)
    # confidence downgrade when GMM disagrees
    e.loc[(e["gmm_label"] != e["burn_type"]) & (e["type_confidence"] == "high"), "type_confidence"] = "medium"
    e["conflict_event"] = e["conflict_frac"] > 0.5
    keep = ["event_id", "start", "end", "year", "doy", "n_det", "duration_days", "extent_km", "frp_mean", "frp_max",
            "frp_sum", "night_frac", "lat", "lon", "oblast", "lc_crop", "lc_grass", "lc_tree", "lc_wetland",
            "lc_built", "lc_other", "recur_frac", "conflict_frac", "peat_frac", "conflict_event", "burn_type",
            "burn_type_wc_only", "type_confidence",
            "gmm_cluster", "gmm_label"]
    ev = e[keep].copy()
    for c in ev.select_dtypes("float64").columns:
        ev[c] = ev[c].astype("float32")
    for c in ("oblast", "burn_type", "type_confidence", "gmm_label"):
        ev[c] = ev[c].astype("category")
    ev["start"] = ev["start"].dt.tz_localize(None)
    ev["end"] = ev["end"].dt.tz_localize(None)
    ev.to_parquet(os.path.join(INTERIM, "events_all.parquet"), index=False, compression="zstd")
    v[["datetime_utc", "lat", "lon", "satellite", "event_id"]].to_parquet(
        os.path.join(INTERIM, "viirs_event_members.parquet"), index=False)
    ct.to_csv(os.path.join(PROC, "events_gmm_vs_rules_crosstab.csv"))
    prof.to_csv(os.path.join(PROC, "events_gmm_cluster_profiles.csv"))
    pd.Series(stats).to_csv(os.path.join(PROC, "events_gmm_agreement.csv"))
    print("events", len(ev), "MB", round(os.path.getsize(os.path.join(INTERIM, "events_all.parquet")) / 1e6, 2))
    print(ev["burn_type"].value_counts().to_string())
    print(ev.groupby("burn_type", observed=True)[["n_det", "duration_days", "frp_mean", "night_frac"]].median().round(2).to_string())
    print("detections by type:", ev.groupby("burn_type", observed=True)["n_det"].sum().to_dict())
    print("GMM vs rules:", stats)
    print(ct.to_string())
    print(prof.round(2).to_string())
    print("size distribution n_det:", ev["n_det"].describe(percentiles=[.5, .9, .99]).round(1).to_dict())
    print(f"t={time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
