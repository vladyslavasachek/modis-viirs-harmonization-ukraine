"""Stage 2: MODIS -> VIIRS S-NPP-equivalent harmonization.

Unit: model cell (0.25 or 0.5 deg) x year x week-of-year (0..51; week 51 absorbs days 358-366) x day/night.
Target: S-NPP detection count (standard processing). Predictors from MODIS (Terra+Aqua, or Terra-only).
Models: B0 global ratio; B1 ratio by month x day/night with empirical-Bayes shrinkage to B0;
M1 Negative-Binomial GLM with explicit MODIS==0 ("VIIRS-only small fires") terms.
Validation: leave-one-year-out (LOYO) on the pre-invasion years 2012-2021.
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INTERIM = os.path.join(ROOT, "data", "interim")
PROC = os.path.join(ROOT, "data", "processed")
LC = ["crop", "grass", "tree", "wetland"]
TRAIN_YEARS = list(range(2012, 2022))      # main fit: pre-invasion
SENS_YEARS = list(range(2012, 2026))       # sensitivity: incl. 2022-2025, conflict-flagged detections excluded
SNPP_START = pd.Timestamp("2012-01-25")


# ----------------------------------------------------------------------------- data
def load_det() -> pd.DataFrame:
    d = pd.read_parquet(os.path.join(INTERIM, "detections_all.parquet"),
                        columns=["satellite", "processing", "date", "daynight", "lat", "lon", "frp",
                                 "in_ukraine", "conflict_zone", "oblast", "cell01", "confidence", "datetime_utc",
                                 "landcover"])
    d = d[d["in_ukraine"]].copy()
    d["year"] = d["date"].dt.year.astype("int16")
    d["woy"] = np.minimum((d["date"].dt.dayofyear - 1) // 7, 51).astype("int8")
    d["dn"] = (d["daynight"] == "N").astype("int8")
    return d


def cell_id(lat, lon, res):
    k = int(round(1 / res * 100))  # cells per 100 deg
    return (np.floor(np.asarray(lat) * k / 100).astype(int) * 10000 + np.floor(np.asarray(lon) * k / 100).astype(int))


def week_month(year, woy):
    mid = pd.to_datetime(year.astype(str) + "-01-01") + pd.to_timedelta(woy.astype(int) * 7 + 3, "D")
    return mid.dt.month.values if hasattr(mid, "dt") else mid.month.values


def outage_weeks(d: pd.DataFrame) -> pd.DataFrame:
    """Weeks with >=1 S-NPP outage day: no S-NPP detection in Ukraine while NOAA-20 >=5 or MODIS >=3."""
    daily = d.groupby(["date", "satellite"]).size().unstack(fill_value=0)
    for c in ("SNPP", "NOAA20", "Terra", "Aqua"):
        if c not in daily:
            daily[c] = 0
    out = (daily["SNPP"] == 0) & ((daily["NOAA20"] >= 5) | (daily["Terra"] + daily["Aqua"] >= 3))
    out |= daily.index < SNPP_START
    od = daily.index[out]
    w = pd.DataFrame({"year": od.year, "woy": np.minimum((od.dayofyear - 1) // 7, 51)}).drop_duplicates()
    # also the partial first S-NPP weeks of 2012
    w = pd.concat([w, pd.DataFrame({"year": 2012, "woy": range(0, 4)})]).drop_duplicates()
    w["outage"] = True
    return w, od


def lc_shares(res) -> pd.DataFrame:
    s = pd.read_parquet(os.path.join(PROC, "cell01_landcover_shares.parquet"))
    lat = (s["cell01"] // 1000) / 10 + 0.05
    lon = (s["cell01"] % 1000) / 10 + 0.05
    s["cell"] = cell_id(lat, lon, res)
    cols = [c for c in LC if c in s]
    return s.groupby("cell")[cols].mean().reset_index()


def build_table(d: pd.DataFrame, res: float, years, modis_sats=("Terra", "Aqua"), exclude_conflict=False,
                target_sat="SNPP", processing=("standard",)):
    """Full grid cell x year x woy x dn with MODIS predictors and target counts."""
    d = d[d["year"].isin(years) & d["processing"].isin(processing)]
    if exclude_conflict:
        d = d[~d["conflict_zone"]]
    d = d.assign(cell=cell_id(d["lat"], d["lon"], res))
    keys = ["cell", "year", "woy", "dn"]
    m = d[d["satellite"].isin(modis_sats)].groupby(keys).agg(m_n=("frp", "size"), m_frp=("frp", "sum"))
    v = d[d["satellite"] == target_sat].groupby(keys).size().rename("v_n")
    return m, v


def full_grid(cells, years, m, v=None):
    idx = pd.MultiIndex.from_product([cells, years, range(52), (0, 1)], names=["cell", "year", "woy", "dn"])
    g = pd.DataFrame(index=idx).join(m, how="left")
    if v is not None:
        g = g.join(v, how="left")
    g = g.fillna(0).reset_index()
    g["month"] = week_month(g["year"], g["woy"])
    return g


def prepare(g, shares, cell_oblast, outage=None):
    g = g.merge(shares, on="cell", how="left").fillna({c: 0 for c in LC})
    g = g.merge(cell_oblast, on="cell", how="left")
    if outage is not None:
        g = g.merge(outage, on=["year", "woy"], how="left")
        g["outage"] = g["outage"].fillna(False).astype(bool)
    return g


# ----------------------------------------------------------------------------- models
def design(g: pd.DataFrame) -> np.ndarray:
    n = len(g)
    dn = g["dn"].values.astype(float)
    mon = g["month"].values
    M = np.zeros((n, 11))
    M[np.arange(n)[mon > 1], mon[mon > 1] - 2] = 1.0
    lm = np.log1p(g["m_n"].values)
    lf = np.log1p(g["m_frp"].values)
    z = (g["m_n"].values == 0).astype(float)
    lc = g[LC].values
    cols = [np.ones(n)[:, None], dn[:, None], M, M * dn[:, None], lc, lm[:, None], (lm * dn)[:, None],
            lf[:, None], z[:, None], (z * dn)[:, None], M * z[:, None], lc * z[:, None]]
    return np.hstack(cols)


class B0:
    name = "B0"

    def fit(self, g):
        self.r = g["v_n"].sum() / g["m_n"].sum()
        return self

    def predict(self, g):
        return self.r * g["m_n"].values


class B1:
    """Ratio by month x day/night, shrunk to the global ratio (prior strength K MODIS detections)."""
    name = "B1"
    K = 200.0

    def fit(self, g):
        r0 = g["v_n"].sum() / g["m_n"].sum()
        s = g.groupby(["month", "dn"])[["v_n", "m_n"]].sum()
        self.r = ((s["v_n"] + self.K * r0) / (s["m_n"] + self.K)).rename("r")
        self.r0 = r0
        return self

    def predict(self, g):
        r = g[["month", "dn"]].merge(self.r.reset_index(), on=["month", "dn"], how="left")["r"].fillna(self.r0)
        return r.values * g["m_n"].values


class M1:
    name = "M1"

    def __init__(self, alpha=None):
        self.alpha = alpha

    start = None  # class-level warm start shared across fits (LOYO / bootstrap)

    @staticmethod
    def irls(X, y, alpha, beta0=None, maxiter=25, tol=1e-7, ridge=1e-6):
        """NB2 GLM (log link, fixed alpha) by IRLS; tiny ridge for numerical stability only."""
        p = X.shape[1]
        beta = np.zeros(p) if beta0 is None else beta0.copy()
        if beta0 is None:
            beta[0] = np.log(max(y.mean(), 1e-3))
        dev_old = np.inf
        for _ in range(maxiter):
            eta = np.clip(X @ beta, -30, 30)
            mu = np.exp(eta)
            w = mu / (1 + alpha * mu)
            z = eta + (y - mu) / mu
            XtW = X.T * w
            beta_new = np.linalg.solve(XtW @ X + ridge * np.eye(p), XtW @ z)
            step = np.abs(beta_new - beta).max()
            beta = beta_new
            if step < 1e-5:
                break
        return beta

    def fit(self, g):
        X = design(g)
        y = g["v_n"].values.astype(float)
        if self.alpha is None:  # moment estimate of NB2 alpha from a Poisson fit (alpha -> 0 = Poisson)
            bp = self.irls(X, y, 1e-8, maxiter=100)
            mu = np.exp(np.clip(X @ bp, -30, 30))
            zz = ((y - mu) ** 2 - y) / mu
            self.alpha = float(max(np.sum(zz * mu) / np.sum(mu ** 2), 1e-3))
            self.params = self.irls(X, y, self.alpha, beta0=bp, maxiter=100)
            M1.start = self.params  # warm start for subsequent LOYO / bootstrap fits
            return self
        self.params = self.irls(X, y, self.alpha, beta0=M1.start, maxiter=100)
        return self

    def predict(self, g, params=None):
        p = self.params if params is None else params
        return np.exp(np.clip(design(g) @ p, -30, 30))


# ----------------------------------------------------------------------------- metrics
LEVELS = {"ukr_month": ["year", "month"], "ukr_year": ["year"], "oblast_month": ["oblast", "year", "month"],
          "cell_month": ["cell", "year", "month"]}


def metrics(df, pred_col="pred", obs_col="v_n"):
    out = {}
    for lev, keys in LEVELS.items():
        a = df.groupby(keys)[[pred_col, obs_col]].sum()
        err = a[pred_col] - a[obs_col]
        out[f"{lev}_WAPE"] = np.abs(err).sum() / a[obs_col].sum()
        out[f"{lev}_bias"] = err.sum() / a[obs_col].sum()
    return out


def loyo(g, years, model_factories):
    """Returns per-model concatenated held-out predictions."""
    preds = {}
    for mk, fac in model_factories.items():
        parts = []
        for y in years:
            tr = g[(g["year"] != y) & g["year"].isin(years)]
            te = g[g["year"] == y].copy()
            m = fac().fit(tr)
            te["pred"] = m.predict(te)
            parts.append(te)
        preds[mk] = pd.concat(parts)
    return preds


# ----------------------------------------------------------------------------- AOI series with bootstrap replicates
SOURCES = {0: "S-NPP observed", 1: "NOAA-20 calibrated", 2: "MODIS-based estimate (Terra+Aqua)",
           3: "MODIS-based estimate (Terra only, lower confidence)"}


def ratio_table(model: "B1") -> np.ndarray:
    """B1 ratios as array [12 months, 2 dn] (missing strata -> global ratio)."""
    R = np.full((12, 2), model.r0)
    for (mo, dn), r in model.r.items():
        R[int(mo) - 1, int(dn)] = r
    return R


PRED_BINS = [0, 1e-9, 3, 10, 30, 100, 300, 1000, 3000, 10000, np.inf]
QS = {"lo95": .025, "lo80": .1, "hi80": .9, "hi95": .975}


SCALES = {"cell": ["cell"], "oblast": ["oblast"], "ukr": []}


def aoi_scale(n_cells: int) -> str:
    """Which calibration to use for an AOI of n 0.5-deg cells (Ukraine = 365, oblasts ~10-30)."""
    return "cell" if n_cells <= 3 else ("oblast" if n_cells <= 80 else "ukr")


def calibrate_by_scale(loyo_pred: pd.DataFrame) -> pd.DataFrame:
    return pd.concat([calibrate_intervals(loyo_pred, {sc: sk}).assign(scale=sc) for sc, sk in SCALES.items()],
                     ignore_index=True)


def calibrate_intervals(loyo_pred: pd.DataFrame, levels=None, min_n=30) -> pd.DataFrame:
    """Empirical (LOYO) predictive intervals for MODIS-based estimates.

    Pools held-out aggregates at Ukraine / oblast / cell x week / month / year, bins them by predicted size,
    and stores quantiles of obs/pred (bin 0: quantiles of obs itself, since pred = 0)."""
    levels = levels or {"ukr": [], "oblast": ["oblast"], "cell": ["cell"]}
    frames = []
    for space, sk in levels.items():
        for time_, tk in (("week", ["year", "woy"]), ("month", ["year", "month"]), ("year", ["year"])):
            a = loyo_pred.groupby(sk + tk)[["pred", "v_n"]].sum()
            frames.append(a.assign(level=f"{space}_{time_}"))
    a = pd.concat(frames)
    a["bin"] = pd.cut(a["pred"], PRED_BINS, right=False, labels=False)
    rows = []
    for b, g in a.groupby("bin"):
        if len(g) < min_n:
            continue   # sparse bin -> filled from neighbouring bin in interval_from_calib
        x = g["v_n"] if b == 0 else g["v_n"] / g["pred"]
        rows.append(dict(bin=int(b), lo=PRED_BINS[int(b)], hi=PRED_BINS[int(b) + 1], n=len(g),
                         **{k: float(np.quantile(x, q)) for k, q in QS.items()}))
    return pd.DataFrame(rows)


def interval_from_calib(P: np.ndarray, calib: pd.DataFrame, scale: str | None = None) -> dict:
    P = np.asarray(P, dtype=float)
    if scale is not None and "scale" in calib:
        calib = calib[calib["scale"] == scale]
    b = np.clip(np.searchsorted(np.array(PRED_BINS), P, side="right") - 1, 0, len(PRED_BINS) - 2)
    c = calib.set_index("bin").reindex(range(len(PRED_BINS) - 1)).ffill().bfill()
    out = {}
    for k in QS:
        q = c[k].values[b]
        out[k] = np.where(b == 0, q, P * q)
    return out


def aoi_components(t: pd.DataFrame, variant: str = "excl", n20_ratio: float = 1.0) -> pd.DataFrame:
    """Per (year, woy, month) AOI totals split by source; t = rows of harmonized_cells for the AOI."""
    v = variant
    t = t.assign(obs=np.where(t["source"] == 0, t[f"snpp_{v}"], 0.0),
                 n20=np.where(t["source"] == 1, t[f"n20_{v}"] * n20_ratio, 0.0),
                 mod_ta=np.where(t["source"] == 2, t[f"h_{v}"], 0.0),
                 mod_t=np.where(t["source"] == 3, t[f"h_{v}"], 0.0),
                 has_ta=(t["source"] == 2), has_t=(t["source"] == 3))
    g = t.groupby(["year", "woy"]).agg(month=("month", "min"), source=("source", "max"), obs=("obs", "sum"),
                                       n20=("n20", "sum"), mod_ta=("mod_ta", "sum"), mod_t=("mod_t", "sum"),
                                       has_ta=("has_ta", "max"), has_t=("has_t", "max"))
    g["point"] = g[["obs", "n20", "mod_ta", "mod_t"]].sum(axis=1)
    return g.reset_index()


def summarize_components(c: pd.DataFrame, keys, calib_ta, calib_t, n20q: dict, scale: str = "ukr") -> pd.DataFrame:
    """Aggregate weekly components to `keys` and attach 80/95% intervals.
    Observed S-NPP: no model interval. NOAA-20 part: ratio bootstrap. MODIS parts: LOYO-calibrated.
    Parts are added quantile-wise (conservative)."""
    a = c.groupby(keys).agg(obs=("obs", "sum"), n20=("n20", "sum"), mod_ta=("mod_ta", "sum"),
                            mod_t=("mod_t", "sum"), point=("point", "sum"), has_ta=("has_ta", "max"),
                            has_t=("has_t", "max"))
    ita = interval_from_calib(a["mod_ta"].values, calib_ta, scale)
    it = interval_from_calib(a["mod_t"].values, calib_t, scale)
    for k in QS:
        a[k] = a["obs"] + a["n20"] * n20q[k] + np.where(a["has_ta"], ita[k], 0) + np.where(a["has_t"], it[k], 0)
    a["share_estimated"] = (a["mod_ta"] + a["mod_t"] + a["n20"]) / a["point"].where(a["point"] > 0)
    return a.reset_index()


def dominant_lc(g: pd.DataFrame) -> np.ndarray:
    s = g[LC].values
    other = 1 - s.sum(axis=1)
    allc = np.column_stack([s, other])
    return np.array(LC + ["other"])[allc.argmax(axis=1)]


class B2(B1):
    """B1 x land-cover: ratio by month x day/night x dominant cell land cover, shrunk to the B1 ratio."""
    name = "B2"
    K2 = 200.0

    def fit(self, g):
        super().fit(g)
        g = g.assign(lc_dom=dominant_lc(g))
        s = g.groupby(["month", "dn", "lc_dom"])[["v_n", "m_n"]].sum().reset_index()
        s = s.merge(self.r.reset_index(), on=["month", "dn"])
        s["r2"] = (s["v_n"] + self.K2 * s["r"]) / (s["m_n"] + self.K2)
        self.r2 = s[["month", "dn", "lc_dom", "r2"]]
        return self

    def predict(self, g):
        k = g[["month", "dn"]].assign(lc_dom=dominant_lc(g))
        r = k.merge(self.r2, on=["month", "dn", "lc_dom"], how="left")["r2"]
        r = r.fillna(pd.Series(super().predict(g) / np.maximum(g["m_n"].values, 1e-9)))
        return r.values * g["m_n"].values
