"""Stage 6: daily upwind fire load for Kyiv and its association with PM2.5 (CAMS via Open-Meteo).

Upwind fire load (UFL) for day t and lag k (k = 0, 1, 2 days):
  UFL_k(t) = sum over detections on day t-k within R_MAX km of Kyiv whose bearing from Kyiv lies within
             +-HALF_SECTOR deg of the direction the wind blows FROM (vector-mean ERA5 100 m wind over days t-k..t)
             of FRP * exp(-distance / L_DECAY).
Fire stream: VIIRS S-NPP (reference sensor) 2012 - 2 Nov 2026, then NOAA-20; Ukraine + any extra FIRMS region
(data/raw/firms_*/, e.g. Belarus). Conflict-affected detections (app default domain) are excluded.
Combined load: UFL = UFL_0 + 0.5 UFL_1 + 0.25 UFL_2.

PM: CAMS European regional ensemble via Open-Meteo Air Quality API (daily mean of hourly PM2.5 at Kyiv).
CAVEAT: CAMS assimilates GFAS fire emissions derived from MODIS/VIIRS FRP, so a CAMS-PM vs FIRMS-FRP association
is partly circular. No ground-station data were used (OpenAQ requires an API key).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, os.path.dirname(__file__))
from harmonize import INTERIM, PROC, ROOT

KYIV = (50.4501, 30.5234)
R_MAX, L_DECAY, HALF_SECTOR = 300.0, 100.0, 30.0
LAG_W = {0: 1.0, 1: 0.5, 2: 0.25}
OUT = os.path.join(PROC, "kyiv")
os.makedirs(OUT, exist_ok=True)


def dist_bearing(lat, lon):
    la0, lo0 = np.radians(KYIV[0]), np.radians(KYIV[1])
    la, lo = np.radians(lat), np.radians(lon)
    dlo = lo - lo0
    a = np.sin((la - la0) / 2) ** 2 + np.cos(la0) * np.cos(la) * np.sin(dlo / 2) ** 2
    d = 6371 * 2 * np.arcsin(np.sqrt(a))
    b = np.degrees(np.arctan2(np.sin(dlo) * np.cos(la), np.cos(la0) * np.sin(la) - np.sin(la0) * np.cos(la) * np.cos(dlo)))
    return d, (b + 360) % 360


def load_weather():
    m = json.load(open(os.path.join(ROOT, "data", "raw", "meteo", "era5_kyiv.json")))
    h = pd.DataFrame(m["hourly"])
    h["time"] = pd.to_datetime(h["time"])
    # wind blowing FROM direction theta: u = -ws sin(theta), v = -ws cos(theta)
    th = np.radians(h["wind_direction_100m"])
    h["u"], h["v"] = -h["wind_speed_100m"] * np.sin(th), -h["wind_speed_100m"] * np.cos(th)
    dly = h.groupby(h["time"].dt.floor("D")).agg(u=("u", "mean"), v=("v", "mean"), ws100=("wind_speed_100m", "mean"),
                                                   ws10=("wind_speed_10m", "mean"), blh=("boundary_layer_height", "mean"))
    d = pd.DataFrame(m["daily"])
    d["time"] = pd.to_datetime(d["time"])
    dly = dly.join(d.set_index("time"))
    dly.index.name = "date"
    return dly


def wind_from(u, v):
    return (np.degrees(np.arctan2(-u, -v)) + 360) % 360


def load_fires():
    det = pd.read_parquet(os.path.join(PROC, "app", "detections.parquet"),
                          columns=["sensor", "datetime_utc", "lat", "lon", "frp", "oblast", "burn_type",
                                   "conflict_flag", "weight"])
    det["lat"], det["lon"], det["frp"] = det["lat"] / 1e4, det["lon"] / 1e4, det["frp"] / 10
    det["datetime_utc"] = pd.to_datetime(det["datetime_utc"].astype("int64") * 60, unit="s")
    snpp_end = det.loc[det["sensor"] == "VIIRS S-NPP", "datetime_utc"].max().floor("D")
    ua = det[((det["sensor"] == "VIIRS S-NPP") |
              ((det["sensor"] == "VIIRS NOAA-20") & (det["datetime_utc"] >= snpp_end + pd.Timedelta(days=1))))
             & ~det["conflict_flag"]]
    ua = pd.DataFrame({"date": ua["datetime_utc"].dt.floor("D"), "lat": ua["lat"].astype(float),
                       "lon": ua["lon"].astype(float), "frp": ua["frp"].astype(float),
                       "region": ua["oblast"].astype(str), "burn_type": ua["burn_type"].astype(str)})
    parts = [ua]
    p = os.path.join(INTERIM, "extra_detections.parquet")
    if os.path.exists(p):
        x = pd.read_parquet(p)
        x = x[(x["satellite"] == "SNPP") | ((x["satellite"] == "NOAA20") & (x["date"] > snpp_end))]
        parts.append(pd.DataFrame({"date": x["date"], "lat": x["lat"], "lon": x["lon"], "frp": x["frp"],
                                   "region": x["region"], "burn_type": x["burn_type"]}))
    f = pd.concat(parts, ignore_index=True)
    f["dist_km"], f["bearing"] = dist_bearing(f["lat"].values, f["lon"].values)
    f = f[f["dist_km"] <= R_MAX].copy()
    f["w"] = f["frp"] * np.exp(-f["dist_km"] / L_DECAY)
    return f


def upwind_load(f, wx):
    rows = []
    days = wx.index
    for k in (0, 1, 2):
        # mean wind vector over days s..s+k (transport window), assigned to target day t = s+k
        u = wx["u"].rolling(k + 1).mean().shift(-k)
        v = wx["v"].rolling(k + 1).mean().shift(-k)
        wd = pd.Series(wind_from(u.values, v.values), index=days, name="wd")
        g = f.merge(wd.rename("wd"), left_on="date", right_index=True, how="inner")
        diff = np.abs((g["bearing"] - g["wd"] + 180) % 360 - 180)
        g = g[diff <= HALF_SECTOR].assign(target=lambda x: x["date"] + pd.Timedelta(days=k), lag=k)
        rows.append(g[["target", "lag", "w", "region", "burn_type", "dist_km"]])
    return pd.concat(rows, ignore_index=True)


def main():
    wx = load_weather()
    f = load_fires()
    print("fires within 300 km:", len(f), f.groupby(f["region"] == "Belarus")["w"].sum().round(0).to_dict())
    up = upwind_load(f, wx)
    up["wl"] = up["w"] * up["lag"].map(LAG_W)
    daily = up.pivot_table(index="target", columns="lag", values="w", aggfunc="sum").reindex(wx.index).fillna(0.0)
    daily.columns = [f"ufl_lag{c}" for c in daily.columns]
    daily["ufl"] = sum(daily[f"ufl_lag{k}"] * w for k, w in LAG_W.items())
    by_type = up.pivot_table(index="target", columns="burn_type", values="wl", aggfunc="sum").reindex(wx.index).fillna(0)
    up["region2"] = np.where(up["region"] == "Belarus", "Belarus", up["region"])
    by_reg = up.pivot_table(index="target", columns="region2", values="wl", aggfunc="sum").reindex(wx.index).fillna(0)

    # PM2.5 daily mean (CAMS)
    a = json.load(open(os.path.join(ROOT, "data", "raw", "aq", "cams_kyiv_hourly.json")))
    pm = pd.DataFrame(a["hourly"])
    pm["time"] = pd.to_datetime(pm["time"])
    pmd = pm.groupby(pm["time"].dt.floor("D")).agg(pm25=("pm2_5", "mean"), pm10=("pm10", "mean"), n=("pm2_5", "count"))
    pmd = pmd[pmd["n"] >= 18].drop(columns="n")
    df = wx.join(daily).join(pmd)
    df["wd100"] = wind_from(df["u"].values, df["v"].values)
    df.index.name = "date"
    df.reset_index().to_parquet(os.path.join(OUT, "kyiv_daily.parquet"), index=False)
    by_type.rename_axis("date").reset_index().to_parquet(os.path.join(OUT, "kyiv_ufl_by_type.parquet"), index=False)
    by_reg.rename_axis("date").reset_index().to_parquet(os.path.join(OUT, "kyiv_ufl_by_region.parquet"), index=False)

    # ---------------- regression: PM anomaly ~ meteorology (+ fire load); chronological validation
    d = df.dropna(subset=["pm25"]).copy()
    d["lpm"] = np.log(d["pm25"])
    # anomaly vs trailing level (days t-365..t-31): removes slow level shifts (CAMS product changes) w/o using the future
    d["level"] = d["lpm"].rolling(335, min_periods=200).median().shift(31)
    d = d.dropna(subset=["level"])
    d["anom"] = d["lpm"] - d["level"]
    d["month"] = d.index.month
    X_met = pd.get_dummies(d["month"], prefix="m", drop_first=True).astype(float)
    X_met["ws100"] = d["ws100"]
    X_met["log_blh"] = np.log(d["blh"].clip(lower=50)).fillna(np.log(d["blh"].median()))
    X_met["precip"] = np.log1p(d["precipitation_sum"])
    X_met["tmax"] = d["temperature_2m_max"]
    X_met["rh"] = d["relative_humidity_2m_mean"]
    X_met["wd_sin"], X_met["wd_cos"] = np.sin(np.radians(d["wd100"])), np.cos(np.radians(d["wd100"]))
    X_fire = X_met.copy()
    for k in (0, 1, 2):
        X_fire[f"log1p_ufl_lag{k}"] = np.log1p(d[f"ufl_lag{k}"])
    y = d["anom"]
    train = d.index < "2022-01-01"
    res = []
    for name, X in (("meteorology only", X_met), ("meteorology + upwind fire load", X_fire)):
        Xc = sm.add_constant(X)
        fit = sm.OLS(y[train], Xc[train]).fit()
        pred = fit.predict(Xc[~train])
        e = y[~train] - pred
        r2 = 1 - (e ** 2).sum() / ((y[~train] - y[~train].mean()) ** 2).sum()
        # skill on high-PM days of the test period
        hi = d.loc[~train, "pm25"] > d["pm25"].quantile(0.9)
        res.append(dict(model=name, train="2013-2021", test="2022-2026", n_train=int(train.sum()),
                        n_test=int((~train).sum()), test_R2=r2, test_RMSE=float(np.sqrt((e ** 2).mean())),
                        test_RMSE_highPM_days=float(np.sqrt((e[hi] ** 2).mean())), n_test_highPM=int(hi.sum())))
    val = pd.DataFrame(res)
    val.to_csv(os.path.join(OUT, "kyiv_regression_validation.csv"), index=False)
    print(val.round(4).to_string(index=False))
    full = sm.OLS(y, sm.add_constant(X_fire)).fit(cov_type="HAC", cov_kwds={"maxlags": 7})
    coef = pd.DataFrame({"coef": full.params, "se_HAC": full.bse, "p": full.pvalues,
                         "ci_lo": full.conf_int()[0], "ci_hi": full.conf_int()[1]})
    coef = coef.loc[[c for c in coef.index if c.startswith("log1p_ufl") or c in ("ws100", "log_blh", "precip", "tmax", "rh")]]
    coef["pct_change_PM_per_doubling_of_(1+UFL)"] = (np.exp(coef["coef"] * np.log(2)) - 1) * 100
    coef.to_csv(os.path.join(OUT, "kyiv_regression_coefficients.csv"))
    print(coef.round(4).to_string())

    # ---------------- high-PM days
    d2 = df.dropna(subset=["pm25"]).copy()
    p90_pm = d2["pm25"].quantile(0.9)
    pos = d2["ufl"] > 0
    p90_ufl = d2.loc[pos, "ufl"].quantile(0.9) if pos.any() else 0
    d2["high_pm"] = d2["pm25"] > p90_pm
    d2["high_ufl"] = d2["ufl"] > p90_ufl
    summ = {"pm25_p90_ugm3": p90_pm, "ufl_p90_of_nonzero_days": p90_ufl,
            "share_all_days_high_ufl": d2["high_ufl"].mean(),
            "share_highPM_days_with_high_ufl": d2.loc[d2["high_pm"], "high_ufl"].mean(),
            "share_highPM_days_with_any_upwind_fire": (d2.loc[d2["high_pm"], "ufl"] > 0).mean(),
            "share_all_days_with_any_upwind_fire": (d2["ufl"] > 0).mean(),
            "n_days": len(d2), "n_high_pm_days": int(d2["high_pm"].sum())}
    for season, months in (("spring (Mar-May)", [3, 4, 5]), ("summer-autumn (Jul-Oct)", [7, 8, 9, 10])):
        s = d2[d2.index.month.isin(months)]
        summ[f"{season}: share highPM days with high ufl"] = s.loc[s["high_pm"], "high_ufl"].mean()
        summ[f"{season}: share all days with high ufl"] = s["high_ufl"].mean()
    summ["lift_highPM_vs_all"] = summ["share_highPM_days_with_high_ufl"] / summ["share_all_days_high_ufl"]
    pd.Series(summ).to_csv(os.path.join(OUT, "kyiv_highpm_summary.csv"))
    print(pd.Series(summ).round(3).to_string())
    sel = d2.index[d2["high_pm"] & d2["high_ufl"]]
    bt = by_type.loc[sel].sum()
    br = by_reg.loc[sel].sum()
    (bt / bt.sum()).sort_values(ascending=False).rename("share_of_upwind_load").to_csv(
        os.path.join(OUT, "kyiv_highpm_load_by_type.csv"))
    (br / br.sum()).sort_values(ascending=False).rename("share_of_upwind_load").to_csv(
        os.path.join(OUT, "kyiv_highpm_load_by_region.csv"))
    print("on high-PM & high-load days, by type:", (bt / bt.sum()).round(3).sort_values(ascending=False).to_dict())
    print("by region:", (br / br.sum()).round(3).sort_values(ascending=False).head(8).to_dict())
    # episodes for case studies: top weeks by combined load with PM above p90
    ep = d2[d2["high_pm"] & d2["high_ufl"]].sort_values("ufl", ascending=False)
    print(ep[["pm25", "ufl", "ufl_lag0", "ufl_lag1", "ws100", "wd100"]].head(15).round(1).to_string())


if __name__ == "__main__":
    main()
