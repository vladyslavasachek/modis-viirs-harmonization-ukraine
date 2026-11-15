"""Stage 6b: Kyiv upwind fire load vs ground-station PM2.5 (OpenAQ, low-cost AirGradient sensors, 2023-09 -> 2026-11).

Station composite: daily median of PM2.5 over all OpenAQ PM2.5 sensors within 25 km of Kyiv with >= 75% hourly
coverage that day. Low-cost optical sensors (no humidity correction applied by us): absolute values are uncertain,
day-to-day variation is informative. Same regression design as kyiv_air.py but without the trailing-level term
(record too short); month fixed effects; chronological split train < 2025-07-01, test >= 2025-07-01.
CAMS (Open-Meteo) is evaluated on the same days as a secondary series.
"""
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, os.path.dirname(__file__))
from harmonize import PROC

KY = os.path.join(PROC, "kyiv")
SPLIT = "2025-07-01"


def design(d, fire):
    X = pd.get_dummies(d.index.month, prefix="m", drop_first=True).astype(float)
    X.index = d.index
    X["ws100"] = d["ws100"]
    X["log_blh"] = np.log(d["blh"].clip(lower=50)).fillna(np.log(d["blh"].median()))
    X["precip"] = np.log1p(d["precipitation_sum"])
    X["tmax"] = d["temperature_2m_max"]
    X["rh"] = d["relative_humidity_2m_mean"]
    X["wd_sin"], X["wd_cos"] = np.sin(np.radians(d["wd100"])), np.cos(np.radians(d["wd100"]))
    if fire:
        for k in (0, 1, 2):
            X[f"log1p_ufl_lag{k}"] = np.log1p(d[f"ufl_lag{k}"])
    return sm.add_constant(X)


def analyse(d, ycol, label):
    d = d.dropna(subset=[ycol]).copy()
    y = np.log(d[ycol].clip(lower=1))
    tr = d.index < SPLIT
    rows = []
    for name, fire in (("meteorology only", False), ("meteorology + upwind fire load", True)):
        X = design(d, fire)
        fit = sm.OLS(y[tr], X[tr]).fit()
        e = y[~tr] - fit.predict(X[~tr])
        hi = d.loc[~tr, ycol] > d[ycol].quantile(0.9)
        rows.append(dict(pm_series=label, model=name, n_train=int(tr.sum()), n_test=int((~tr).sum()),
                         test_R2=1 - (e ** 2).sum() / ((y[~tr] - y[~tr].mean()) ** 2).sum(),
                         test_RMSE=float(np.sqrt((e ** 2).mean())),
                         test_RMSE_highPM=float(np.sqrt((e[hi] ** 2).mean())), n_test_highPM=int(hi.sum())))
    full = sm.OLS(y, design(d, True)).fit(cov_type="HAC", cov_kwds={"maxlags": 7})
    cf = pd.DataFrame({"pm_series": label, "coef": full.params, "p": full.pvalues,
                       "ci_lo": full.conf_int()[0], "ci_hi": full.conf_int()[1]}).loc[
        [f"log1p_ufl_lag{k}" for k in (0, 1, 2)]]
    cf["pct_PM_per_doubling"] = (2 ** cf["coef"] - 1) * 100
    p90, u90 = d[ycol].quantile(0.9), d.loc[d["ufl"] > 0, "ufl"].quantile(0.9)
    hp, hu = d[ycol] > p90, d["ufl"] > u90
    hs = dict(pm_series=label, n_days=len(d), p90=p90, share_highPM_with_highUFL=hu[hp].mean(),
              share_all_highUFL=hu.mean(), lift=hu[hp].mean() / hu.mean() if hu.mean() > 0 else np.nan)
    return pd.DataFrame(rows), cf, hs


def main():
    k = pd.read_parquet(os.path.join(KY, "kyiv_daily.parquet")).set_index("date")
    oq = pd.read_csv(os.path.join(KY, "openaq_kyiv_daily.csv"), parse_dates=["date"])
    oq = oq[(oq["param"] == "pm25") & (oq["coverage_pct"] >= 75) & (oq["value"] > 0) & (oq["value"] < 500)]
    st = oq.groupby("date").agg(pm25_station=("value", "median"), n_stations=("value", "size"))
    d = k.join(st, how="inner")
    d = d[d.index >= "2023-09-11"]
    print(f"station days {len(d)}, {d.index.min():%Y-%m-%d}..{d.index.max():%Y-%m-%d}, "
          f"median stations/day {d['n_stations'].median():.0f}")
    print("Spearman station vs CAMS PM2.5:", round(d[["pm25_station", "pm25"]].corr("spearman").iloc[0, 1], 3),
          "| station vs UFL:", round(d[["pm25_station", "ufl"]].corr("spearman").iloc[0, 1], 3),
          "| CAMS vs UFL (same days):", round(d[["pm25", "ufl"]].corr("spearman").iloc[0, 1], 3))
    vals, cfs, hss = [], [], []
    for col, label in (("pm25_station", "OpenAQ stations (low-cost)"), ("pm25", "CAMS model, same days")):
        v, c, h = analyse(d, col, label)
        vals.append(v), cfs.append(c), hss.append(h)
    V, C, Hs = pd.concat(vals), pd.concat(cfs), pd.DataFrame(hss)
    V.to_csv(os.path.join(KY, "kyiv_station_regression_validation.csv"), index=False)
    C.to_csv(os.path.join(KY, "kyiv_station_regression_coefficients.csv"))
    Hs.to_csv(os.path.join(KY, "kyiv_station_highpm_summary.csv"), index=False)
    d[["pm25_station", "n_stations", "pm25", "ufl", "ufl_lag0", "ufl_lag1", "ufl_lag2"]].reset_index().to_parquet(
        os.path.join(KY, "kyiv_station_daily.parquet"), index=False)
    pd.set_option("display.width", 220)
    print(V.round(3).to_string(index=False)); print(C.round(4).to_string()); print(Hs.round(3).to_string(index=False))
    print("\nTop station PM2.5 days:")
    print(d.sort_values("pm25_station", ascending=False)[["pm25_station", "n_stations", "pm25", "ufl", "ufl_lag0",
                                                           "ufl_lag1", "wd100"]].head(12).round(1).to_string())
    print("\nSep 2024:")
    print(d.loc["2024-08-28":"2024-09-15", ["pm25_station", "pm25", "ufl", "ufl_lag0", "wd100"]].round(1).to_string())


if __name__ == "__main__":
    main()
