"""Stage 2 final (model chosen by LOYO validation = B1: shrunk ratio by month x day/night, 0.5 deg).

- Fit B1 Terra+Aqua and Terra-only on 2012-2021 (pre-invasion), S-NPP outage weeks excluded.
- Year-bootstrap (B=200) -> ratio replicates -> 80/95% intervals; LOYO coverage check of those intervals.
- Sensitivity: fit on 2012-2025 with conflict-flagged detections excluded; compare ratios; out-of-sample
  check of the main fit on 2022-2025.
- NOAA-20 -> S-NPP ratio replicates (year bootstrap).
- Harmonized cell x week x day/night table 2000-2026, Ukraine series, key figure.
"""
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from harmonize import (B0, B1, LC, PROC, QS, ROOT, SENS_YEARS, TRAIN_YEARS, aoi_components, build_table,
                       calibrate_by_scale, cell_id, full_grid, interval_from_calib, lc_shares, load_det, loyo,
                       metrics, outage_weeks, prepare, ratio_table, summarize_components)
from stage2_validate import cell_oblast_map
from etl import CONFLICT_ZONES

# "excl" variant = consistent spatial domain: open-ended conflict oblasts removed in ALL years (otherwise the
# series shows an artificial 2022 drop); time-limited zones (north, Feb-Apr 2022) removed only in their window.
OPEN_ENDED = [o for o, (a, b) in CONFLICT_ZONES.items() if b is None]

RES = 0.5
NBOOT = 200
FIG = os.path.join(ROOT, "figures")
AQUA_START_YW = 2002 * 100 + 26   # 2002-07-04 falls in week-of-year 26 -> Terra-only through this week
rng = np.random.default_rng(42)


def table(d, sats, years, outage, cells, shares, cob, exclude_conflict=False, target="SNPP"):
    m, v = build_table(d, RES, years, modis_sats=sats, exclude_conflict=exclude_conflict, target_sat=target)
    g = prepare(full_grid(cells, years, m, v), shares, cob, outage)
    return g[~g["outage"]].copy()


def boot_ratios(g, years, B):
    """[B+1, 12, 2] ratio tables; index 0 = fit on all years. Resampling unit = year."""
    by = {y: g[g["year"] == y] for y in years}
    R = [ratio_table(B1().fit(g))]
    for _ in range(B):
        ys = rng.choice(years, size=len(years), replace=True)
        R.append(ratio_table(B1().fit(pd.concat([by[y] for y in ys]))))
    return np.stack(R)


def coverage(g, years, B=200):
    rows = []
    for y in years:
        tr_years = [x for x in years if x != y]
        R = boot_ratios(g[g["year"].isin(tr_years)], tr_years, B)
        te = g[g["year"] == y]
        strat = (te["month"].values - 1) * 2 + te["dn"].values
        Rf = R.reshape(R.shape[0], 24)
        pred = te["m_n"].values[:, None] * Rf[1:, strat].T  # rows x B
        df = pd.concat([te[["year", "month"]].reset_index(drop=True),
                        pd.DataFrame(pred, columns=[f"r{i}" for i in range(B)])], axis=1)
        df["v_n"] = te["v_n"].values
        for lev, keys in (("ukr_month", ["year", "month"]), ("ukr_year", ["year"]), ("oblast_month", None)):
            if keys is None:
                df2 = df.assign(oblast=te["oblast"].values)
                a = df2.groupby(["oblast", "year", "month"]).sum(numeric_only=True)
            else:
                a = df.groupby(keys).sum(numeric_only=True)
            rr = a[[f"r{i}" for i in range(B)]].values
            for lvl, (lo, hi) in (("80", (.1, .9)), ("95", (.025, .975))):
                inside = (a["v_n"].values >= np.quantile(rr, lo, axis=1)) & (a["v_n"].values <= np.quantile(rr, hi, axis=1))
                rows.append({"level": lev, "interval": lvl, "n": len(a), "covered": int(inside.sum())})
    c = pd.DataFrame(rows).groupby(["level", "interval"])[["n", "covered"]].sum()
    c["coverage"] = c["covered"] / c["n"]
    return c.reset_index()


def main():
    t0 = time.time()
    d = load_det()
    d["conflict_zone"] = d["conflict_zone"] | d["oblast"].isin(OPEN_ENDED)   # domain-consistent exclusion
    outage, _ = outage_weeks(d)
    cells = np.unique(cell_id(d["lat"], d["lon"], RES))
    shares, cob = lc_shares(RES), cell_oblast_map(d, RES)

    # ---------------- main fits + bootstrap
    R, G = {}, {}
    for key, sats in (("TA", ("Terra", "Aqua")), ("T", ("Terra",))):
        G[key] = table(d, sats, TRAIN_YEARS, outage, cells, shares, cob)
        R[key] = boot_ratios(G[key], TRAIN_YEARS, NBOOT)
    print(f"fits done t={time.time() - t0:.0f}s", flush=True)

    # ---------------- parameter-only (year-bootstrap) interval coverage: diagnostic, shows it is too narrow
    if os.environ.get("RUN_PARAM_COVERAGE", "0") == "1":   # slow (~2.5 min); result is committed
        cov = pd.concat([coverage(G[k], TRAIN_YEARS).assign(model=f"B1_{k}") for k in ("TA", "T")])
        cov.to_csv(os.path.join(PROC, "validation_coverage_param_bootstrap.csv"), index=False)
        print("param-bootstrap coverage:\n", cov.round(3).to_string(index=False), flush=True)

    # ---------------- predictive intervals calibrated on LOYO residuals + cross-validated coverage
    calib, covrows = {}, []
    for k in ("TA", "T"):
        lp = loyo(G[k], TRAIN_YEARS, {"B1": B1})["B1"]
        calib[k] = calibrate_by_scale(lp)
        calib[k].to_csv(os.path.join(PROC, f"interval_calibration_B1_{k}.csv"), index=False)
        for y in TRAIN_YEARS:   # calibrate without year y, evaluate on year y
            cy = calibrate_by_scale(lp[lp["year"] != y])
            te = lp[lp["year"] == y]
            for space, sk in (("ukr", []), ("oblast", ["oblast"]), ("cell", ["cell"])):
                for tm, tk in (("week", ["year", "woy"]), ("month", ["year", "month"]), ("year", ["year"])):
                    a = te.groupby(sk + tk)[["pred", "v_n"]].sum()
                    a = a[(a["pred"] > 0) | (a["v_n"] > 0)]
                    iv = interval_from_calib(a["pred"].values, cy, space)
                    for lvl, (lo, hi) in (("80", ("lo80", "hi80")), ("95", ("lo95", "hi95"))):
                        ins = (a["v_n"].values >= iv[lo]) & (a["v_n"].values <= iv[hi])
                        covrows.append(dict(model=f"B1_{k}", level=f"{space}_{tm}", interval=lvl, n=len(a),
                                            covered=int(ins.sum())))
    cv = pd.DataFrame(covrows).groupby(["model", "level", "interval"])[["n", "covered"]].sum().reset_index()
    cv["coverage"] = cv["covered"] / cv["n"]
    cv.to_csv(os.path.join(PROC, "validation_interval_coverage.csv"), index=False)
    print(cv.pivot_table(index=["model", "level"], columns="interval", values="coverage").round(3).to_string(), flush=True)

    # ---------------- sensitivity: 2012-2025, conflict-flagged detections excluded
    gs = table(d, ("Terra", "Aqua"), SENS_YEARS, outage, cells, shares, cob, exclude_conflict=True)
    Rs = ratio_table(B1().fit(gs))
    sens = pd.DataFrame([{"month": mo + 1, "night": dn, "ratio_main_2012_2021": R["TA"][0, mo, dn],
                          "boot95_lo": np.quantile(R["TA"][1:, mo, dn], .025),
                          "boot95_hi": np.quantile(R["TA"][1:, mo, dn], .975),
                          "ratio_sens_2012_2025_excl": Rs[mo, dn]} for mo in range(12) for dn in (0, 1)])
    sens["sens_inside_main_95"] = sens["ratio_sens_2012_2025_excl"].between(sens["boot95_lo"], sens["boot95_hi"])
    sens["rel_change"] = sens["ratio_sens_2012_2025_excl"] / sens["ratio_main_2012_2021"] - 1
    sens.to_csv(os.path.join(PROC, "b1_ratio_sensitivity.csv"), index=False)
    print(sens.round(3).to_string(index=False), flush=True)
    oos = gs[gs["year"] >= 2022].copy()
    strat = (oos["month"].values - 1) * 2 + oos["dn"].values
    oos["pred"] = oos["m_n"].values * R["TA"][0].reshape(24)[strat]
    srows = [dict(metrics(oos), model="B1_main(2012-21)_out_of_sample_2022-25_excl")]
    for k, p in loyo(gs, SENS_YEARS, {"B0": B0, "B1": B1}).items():
        srows.append(dict(metrics(p), model=f"{k}_LOYO_2012-25_excl"))
    st = pd.DataFrame(srows)
    st.to_csv(os.path.join(PROC, "validation_sensitivity_2022_2025.csv"), index=False)
    print(st[["model", "ukr_month_WAPE", "ukr_month_bias", "ukr_year_WAPE", "ukr_year_bias", "oblast_month_WAPE"]]
          .round(3).to_string(index=False), flush=True)

    # ---------------- NOAA-20 -> S-NPP ratio replicates (year bootstrap, 2018-2025)
    n20_years = list(range(2018, 2026))
    gn = table(d, ("NOAA20",), n20_years, outage, cells, shares, cob)
    gn = gn[~((gn["year"] == 2018) & (gn["woy"] < 13))]
    yr = gn.groupby("year")[["v_n", "m_n"]].sum()
    n20_ratio = yr["v_n"].sum() / yr["m_n"].sum()
    n20_boot = np.array([(lambda s: s["v_n"] / s["m_n"])(yr.loc[rng.choice(n20_years, len(n20_years))].sum())
                         for _ in range(NBOOT)])

    # ---------------- harmonized cell x week x dn table, 2000-2026
    years = list(range(2000, 2027))
    dd = d.assign(cell=cell_id(d["lat"], d["lon"], RES))
    keys = ["cell", "year", "woy", "dn"]
    parts = {}
    ta, te_ = dd["satellite"].isin(["Terra", "Aqua"]), dd["satellite"] == "Terra"
    for nm, mask in (("m", ta), ("mT", te_), ("mx", ta & ~dd["conflict_zone"]), ("mTx", te_ & ~dd["conflict_zone"])):
        a = dd[mask].groupby(keys).agg(n=("frp", "size"), frp=("frp", "sum"))
        parts[f"{nm}_n"], parts[f"{nm}_frp"] = a["n"], a["frp"]
    for sat, nm in (("SNPP", "snpp"), ("NOAA20", "n20")):
        s = dd[dd["satellite"] == sat]
        parts[f"{nm}_all"] = s.groupby(keys).size()
        parts[f"{nm}_excl"] = s[~s["conflict_zone"]].groupby(keys).size()
    t = prepare(full_grid(cells, years, pd.DataFrame(parts)), shares, cob, outage)
    last = d["date"].max()
    yw = t["year"] * 100 + t["woy"]
    t = t[(yw >= 2000 * 100 + 43) & (yw <= last.year * 100 + min((last.dayofyear - 1) // 7, 51))].copy()
    snpp_last = d.loc[d["satellite"] == "SNPP", "date"].max()
    slw = snpp_last.year * 100 + min((snpp_last.dayofyear - 1) // 7, 51)
    yw = t["year"] * 100 + t["woy"]
    n20_avail = yw >= 2018 * 100 + 13
    snpp_era = (yw >= 2012 * 100 + 3) & (yw < slw)
    t["source"] = np.select(
        [snpp_era & ~t["outage"], (snpp_era & t["outage"] & n20_avail) | (yw >= slw), yw <= AQUA_START_YW],
        [0, 1, 3], default=2).astype("int8")
    strat = (t["month"].values - 1) * 2 + t["dn"].values
    pred_ta = t["m_n"].values * R["TA"][0].reshape(24)[strat]
    pred_t = t["mT_n"].values * R["T"][0].reshape(24)[strat]
    t["pred_modis"] = pred_ta.astype("float32")   # MODIS-based estimate for every week (incl. overlap) for comparison
    pred_tax = t["mx_n"].values * R["TA"][0].reshape(24)[strat]
    pred_tx = t["mTx_n"].values * R["T"][0].reshape(24)[strat]
    for v, pa, pt in (("all", pred_ta, pred_t), ("excl", pred_tax, pred_tx)):
        t[f"h_{v}"] = np.select([t["source"] == 0, t["source"] == 1, t["source"] == 2, t["source"] == 3],
                                [t[f"snpp_{v}"], t[f"n20_{v}"] * n20_ratio, pa, pt]).astype("float32")
    t = t[(t[["m_n", "snpp_all", "n20_all"]].sum(axis=1) > 0)]   # all-zero rows contribute 0 to every sum
    for c in ("m_n", "mT_n", "mx_n", "mTx_n", "snpp_all", "snpp_excl", "n20_all", "n20_excl"):
        t[c] = t[c].astype("int16")
    t = t.drop(columns=["outage", "m_frp", "mT_frp", "mx_frp", "mTx_frp"] + LC)
    t["month"] = t["month"].astype("int8")
    t.to_parquet(os.path.join(PROC, "harmonized_cells.parquet"), index=False, compression="zstd", compression_level=19)
    n20q = {k: float(np.quantile(n20_boot, q) / n20_ratio) for k, q in QS.items()}
    np.savez_compressed(os.path.join(PROC, "harmonization_params.npz"), R_ta=R["TA"], R_t=R["T"],
                        n20_boot=n20_boot, n20_ratio=n20_ratio, res=RES)
    pd.Series(n20q).to_csv(os.path.join(PROC, "noaa20_ratio_rel_quantiles.csv"))
    print("harmonized rows", len(t), "MB", round(os.path.getsize(os.path.join(PROC, "harmonized_cells.parquet")) / 1e6, 2),
          "n20 ratio", round(n20_ratio, 4), f"t={time.time() - t0:.0f}s", flush=True)

    # ---------------- Ukraine series (conflict-flagged excluded from 2022) + key figure
    anns = {}
    for v in ("all", "excl"):
        comp = aoi_components(t, v, n20_ratio)
        comp.to_parquet(os.path.join(PROC, f"ukraine_week_components_{v}.parquet"), index=False)
        summarize_components(comp, ["year", "woy"], calib["TA"], calib["T"], n20q).to_parquet(
            os.path.join(PROC, f"ukraine_week_harmonized_{v}.parquet"), index=False)
        anns[v] = summarize_components(comp, ["year"], calib["TA"], calib["T"], n20q).set_index("year")
    ann = anns["all"].join(d.groupby(["year", "satellite"]).size().unstack(fill_value=0))
    ann["point_excl"] = anns["excl"]["point"]
    ann["modis_estimate"] = t.groupby("year")["pred_modis"].sum()
    ann.to_csv(os.path.join(PROC, "ukraine_annual_harmonized.csv"))
    print(ann.drop(columns=["obs", "n20", "mod_ta", "mod_t"]).round(2).to_string(), flush=True)

    a = ann.loc[2001:2025]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.plot(a.index, a["Terra"] + a["Aqua"], color="#555555", marker="o", ms=3, label="Raw MODIS detections (Terra+Aqua)")
    s = a.loc[2012:]
    ax.plot(s.index, s["SNPP"], color="#d95f02", marker="o", ms=3, lw=1, label="Raw VIIRS S-NPP detections")
    mdl = a.loc[:2011]
    ax.fill_between(mdl.index, mdl["lo95"], mdl["hi95"], color="#1f78b4", alpha=0.15, label="95% interval")
    ax.fill_between(mdl.index, mdl["lo80"], mdl["hi80"], color="#1f78b4", alpha=0.3, label="80% interval")
    ax.plot(a.index, a["point"], color="#1f78b4", lw=2.2, label="Harmonized, S-NPP-equivalent")
    lp = pd.read_parquet(os.path.join(PROC, "loyo_ukr_week.parquet"))
    lpa = lp[lp["model"] == "B1"].groupby("year")["pred"].sum()
    ax.plot(lpa.index, lpa.values, color="#1f78b4", ls=":", lw=1.5,
            label="B1 leave-one-year-out estimate from MODIS (2012–21 check)")
    ax.axvspan(2000.5, 2002.5, color="#eeeeee", zorder=0)
    ymax = ax.get_ylim()[1]
    ax.text(2001.5, ymax * 0.98, "Terra-only\nlower confidence", ha="center", va="top", fontsize=8)
    ax.axvline(2011.5, color="k", lw=0.6, ls="--")
    ax.text(2011.6, ymax * 0.98, "S-NPP observations from 2012", fontsize=8, va="top")
    ax.set_ylabel("hot-spot detections / year (S-NPP-equivalent)")
    ax.set_title("Ukraine: the 2012 sensor jump removed by harmonization\n"
                 "2001–2011: MODIS-based estimate, LOYO-calibrated 80/95% intervals · 2012+: S-NPP observed "
                 "(outage weeks filled with calibrated NOAA-20) · "
                 "all of Ukraine (conflict-zone detections included)", fontsize=10)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "02_harmonized_vs_raw_ukraine.png"), dpi=150)
    print(f"done t={time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
