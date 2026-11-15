"""Stage 1 aggregates, raw annual plot, MODIS version step check."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INTERIM = os.path.join(ROOT, "data", "interim")
PROC = os.path.join(ROOT, "data", "processed")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)


def main():
    det = pd.read_parquet(os.path.join(INTERIM, "detections_all.parquet"))
    ukr = det[det["in_ukraine"]].copy()

    # Base aggregate: 0.1 deg cell x UTC day x day/night x sensor_product x processing x conflict flag
    keys = ["cell01", "date", "daynight", "sensor_product", "satellite", "processing", "conflict_zone"]
    for lc in ("crop", "grass", "tree", "wetland"):
        ukr[f"n_{lc}"] = (ukr["landcover"] == lc).astype("int32")
    agg = ukr.groupby(keys, observed=True).agg(
        n=("frp", "size"), frp_sum=("frp", "sum"), frp_max=("frp", "max"),
        n_crop=("n_crop", "sum"), n_grass=("n_grass", "sum"),
        n_tree=("n_tree", "sum"), n_wetland=("n_wetland", "sum"),
    ).reset_index()
    for c in ("n", "n_crop", "n_grass", "n_tree", "n_wetland"):
        agg[c] = agg[c].astype("int32")
    for c in ("frp_sum", "frp_max"):
        agg[c] = agg[c].astype("float32")
    cells = (ukr.groupby(["cell01", "oblast"]).size().reset_index(name="k")
             .sort_values("k").drop_duplicates("cell01", keep="last").set_index("cell01")["oblast"].rename("oblast_mode"))
    agg = agg.merge(cells, on="cell01", how="left")
    for c in ("daynight", "sensor_product", "satellite", "processing", "oblast_mode"):
        agg[c] = agg[c].astype("category")
    # compact dtypes (sensor_product is implied by satellite)
    agg = agg.drop(columns=["sensor_product"])
    agg["cell01"] = agg["cell01"].astype("int32")
    for c in ("n", "n_crop", "n_grass", "n_tree", "n_wetland"):
        agg[c] = agg[c].astype("int16")
    agg["frp_sum"] = agg["frp_sum"].round(1).astype("float32")
    agg["frp_max"] = agg["frp_max"].round(1).astype("float32")
    agg = agg.sort_values(["cell01", "date"])
    agg.to_parquet(os.path.join(PROC, "agg_cell_day.parquet"), index=False, compression="zstd", compression_level=19)

    # Raw annual counts per satellite (standard + nrt, Ukraine, vegetation, conf-filtered)
    ukr["year"] = ukr["date"].dt.year
    ann = ukr.groupby(["year", "sensor_product", "satellite"]).size().rename("n").reset_index()
    ann.to_csv(os.path.join(PROC, "qc_raw_annual_counts.csv"), index=False)
    fig, ax = plt.subplots(figsize=(10, 5))
    styles = {"Terra": ("#888888", "--"), "Aqua": ("#444444", "-"), "SNPP": ("#d95f02", "-"),
              "NOAA20": ("#1b9e77", "-"), "NOAA21": ("#7570b3", "-")}
    for sat, (col, ls) in styles.items():
        s = ann[ann["satellite"] == sat]
        ax.plot(s["year"], s["n"], color=col, ls=ls, marker="o", ms=3, label=f"{sat} ({s['n'].sum():,})")
    ax.axvline(2026.5, color="k", lw=0.5)
    ax.text(2026.4, ax.get_ylim()[1] * 0.95, "2026 partial", ha="right", fontsize=8)
    ax.set_title("Raw annual hot-spot DETECTIONS per satellite, Ukraine (type 0, confidence-filtered)\n"
                 "Detections are observations, not fires; sensors are not directly comparable")
    ax.set_ylabel("detections / year")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "01_raw_annual_counts.png"), dpi=150)

    # MODIS 6.03 -> 61.03 step check: Aqua detections vs S-NPP (same afternoon orbit), standard only
    std = ukr[ukr["processing"] == "standard"]
    a = std[std["satellite"] == "Aqua"].groupby("year").agg(aqua_n=("frp", "size"), aqua_frp_med=("frp", "median"),
                                                             aqua_ver=("version", lambda s: s.mode().iat[0]))
    a["aqua_conf_med"] = std[std["satellite"] == "Aqua"].assign(c=lambda d: pd.to_numeric(d["confidence"])).groupby("year")["c"].median()
    v = std[std["satellite"] == "SNPP"].groupby("year").size().rename("snpp_n")
    chk = a.join(v).loc[2013:2025]
    chk["ratio_snpp_per_aqua"] = chk["snpp_n"] / chk["aqua_n"]
    chk.to_csv(os.path.join(PROC, "qc_modis_version_check.csv"))
    print(chk.round(2).to_string())
    pre = chk.loc[2017:2022, "ratio_snpp_per_aqua"]
    post = chk.loc[2023:2025, "ratio_snpp_per_aqua"]
    print(f"SNPP/Aqua ratio 2017-22 (v6.03): mean {pre.mean():.2f} range {pre.min():.2f}-{pre.max():.2f}; "
          f"2023-25 (v61.03): mean {post.mean():.2f} range {post.min():.2f}-{post.max():.2f}")

    print("agg rows", len(agg), "size MB", os.path.getsize(os.path.join(PROC, "agg_cell_day.parquet")) / 1e6)
    print(ann.pivot(index="year", columns="satellite", values="n").fillna(0).astype(int).to_string())


if __name__ == "__main__":
    main()
