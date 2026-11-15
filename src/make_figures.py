"""Stage 8: slide figures (PNG) from data/processed only. Colour-blind-safe palettes (viridis/cividis, RdBu, Okabe-Ito)."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "app"))
sys.path.insert(0, os.path.join(ROOT, "src"))
import logging
logging.disable(logging.CRITICAL)
import calendar_stats as cs  # noqa: E402
import lib  # noqa: E402

FIG = os.path.join(ROOT, "figures")
PROC = os.path.join(ROOT, "data", "processed")
UNIT = "VIIRS-equivalent hot spots"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), dpi=150)
    plt.close(fig)
    print("saved", name)


def validation_table():
    v = pd.concat([pd.read_csv(os.path.join(PROC, "validation_loyo.csv")),
                   pd.read_csv(os.path.join(PROC, "validation_loyo_extra.csv"))])
    v = v[(v["grid"] == 0.5) & (v["modis"] == "Terra+Aqua") & v["model"].isin(["B0", "B1", "B2", "M1"])]
    v = v.set_index("model").loc[["B0", "B1", "B2", "M1"]]
    cols = {"ukr_month_WAPE": "Ukraine\n× month", "ukr_year_WAPE": "Ukraine\n× year",
            "oblast_month_WAPE": "Oblast\n× month", "cell_month_WAPE": "Cell (0.5°)\n× month"}
    M = v[list(cols)].values * 100
    names = {"B0": "B0 global ratio", "B1": "B1 month × day/night ratio (chosen)",
             "B2": "B2 = B1 × land cover", "M1": "M1 Negative-Binomial GLM"}
    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(M, cmap="viridis_r", aspect="auto", vmin=0, vmax=40)
    ax.set_xticks(range(4), cols.values())
    ax.set_yticks(range(4), [names[m] for m in v.index])
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{M[i, j]:.1f}%", ha="center", va="center",
                    color="white" if M[i, j] > 22 else "black", fontweight="bold" if i == 1 else None)
    ax.set_title("Leave-one-year-out validation 2012–2021\nWAPE of predicting S-NPP from MODIS (lower is better)")
    fig.colorbar(im, ax=ax, label="WAPE (%)")
    save(fig, "03_validation_table.png")


def ukraine_weekly():
    w = lib.aoi_weekly(tuple(sorted(map(int, lib.aoi_cells("Ukraine")))), "excl")
    return w


def calendar_heatmap(w):
    hm = w[w["year"].between(2001, 2025)].pivot_table(index="year", columns="woy", values="point")
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(hm.values, aspect="auto", cmap="viridis", extent=[0.5, 52.5, 2025.5, 2000.5])
    mt = [pd.Timestamp(2001, m, 1).dayofyear / 7 + 0.5 for m in range(1, 13)]
    ax.set_xticks(mt, [pd.Timestamp(2001, m, 1).strftime("%b") for m in range(1, 13)])
    ax.set_ylabel("year")
    ax.set_title("Burning-activity calendar, Ukraine without conflict-affected oblasts (all years)\n"
                 "weekly harmonized activity; 2001–2011 MODIS-based estimates, 2012+ VIIRS S-NPP observed")
    fig.colorbar(im, ax=ax, label=f"{UNIT} / week")
    save(fig, "04_calendar_heatmap_ukraine.png")


def critical_weeks(w):
    cp, thr = cs.critical_periods(w)
    fig, ax = plt.subplots(figsize=(10, 3.8))
    col = np.where(cp["critical"], "#b2182b", "#92c5de")
    ax.bar(cp["woy"] + 1, cp["prob"] * 100, color=col,
           yerr=[(cp["prob"] - cp["lo95"]) * 100, (cp["hi95"] - cp["prob"]) * 100], ecolor="#555", capsize=1.5)
    ax.axhline(50, ls="--", color="#333", lw=1)
    mt = [pd.Timestamp(2001, m, 1).dayofyear / 7 + 0.5 for m in range(1, 13)]
    ax.set_xticks(mt, [pd.Timestamp(2001, m, 1).strftime("%b") for m in range(1, 13)])
    ax.set_ylabel("probability (%)")
    ax.set_title(f"Critical weeks (red): probability that weekly activity exceeds the baseline 75th percentile "
                 f"(≈{thr:,.0f}/week)\nUkraine without conflict-affected oblasts, 2001–2025, 95% bootstrap intervals")
    save(fig, "05_critical_weeks_ukraine.png")


def burn_type_calendar():
    det = lib.aoi_detections(dict(kind="preset", cells=lib.aoi_cells("Ukraine")), "excl")
    det = det[det["weight"] > 0]
    det = det.assign(m=det["datetime_utc"].dt.month, y=det["datetime_utc"].dt.year, bt=det["burn_type"].astype(str))
    a = det[det["y"].between(2012, 2025)].groupby(["m", "bt"])["weight"].sum().unstack(fill_value=0) / 14
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1.2, 1]})
    bottom = np.zeros(12)
    for t in lib.TYPE_ORDER:
        if t in a:
            axes[0].bar(a.index, a[t], bottom=bottom, color=lib.TYPE_COLORS[t], label=t)
            bottom += a[t].values
    axes[0].set_xticks(range(1, 13), [pd.Timestamp(2001, m, 1).strftime("%b") for m in range(1, 13)])
    axes[0].set_ylabel(f"{UNIT} / month (mean 2012–2025)")
    axes[0].set_title("Burn-type calendar (event-based types, VIIRS 2012–2025)")
    axes[0].legend(fontsize=8, frameon=False)
    y = det.groupby(["y", "bt"])["weight"].sum().unstack(fill_value=0)
    y = y.div(y.sum(axis=1), axis=0)
    bottom = np.zeros(len(y))
    for t in lib.TYPE_ORDER:
        if t in y:
            axes[1].bar(y.index, y[t], bottom=bottom, color=lib.TYPE_COLORS[t], width=0.9)
            bottom += y[t].values
    axes[1].axvline(2011.5, color="k", lw=0.8, ls="--")
    axes[1].text(2006, 1.02, "land-cover-only typing\n(lower confidence)", ha="center", fontsize=8)
    axes[1].set_ylim(0, 1.12)
    axes[1].set_title("Type shares by year")
    axes[1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    fig.text(0.5, 0.01, "Labels are model output, not ground truth (manual check: 14/18 agree, 78%, 95% CI 55–91%)",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(os.path.join(FIG, "06_burn_type_calendar.png"), dpi=150)
    plt.close(fig)
    print("saved 06")


def kyiv_case():
    d = pd.read_parquet(os.path.join(PROC, "kyiv", "kyiv_daily.parquet"))
    d["date"] = pd.to_datetime(d["date"])
    c = d[(d["date"] >= "2020-03-28") & (d["date"] <= "2020-04-30")]
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.bar(c["date"], c["ufl"], color="#E69F00", label="Upwind fire load (FRP, distance-weighted, MW)")
    ax.set_ylabel("upwind fire load (MW)")
    ax2 = ax.twinx()
    ax2.plot(c["date"], c["pm25"], color="#0072B2", lw=2.2, label="PM2.5, CAMS model (µg/m³)")
    ax2.set_ylabel("PM2.5 (µg/m³)")
    ax2.spines["right"].set_visible(True)
    ax.axvspan(pd.Timestamp("2020-04-16"), pd.Timestamp("2020-04-18"), color="#999", alpha=0.25)
    ax.set_ylim(0, 5200)
    ax.text(pd.Timestamp("2020-04-17"), 4400, "Kyiv: world's most\npolluted city\n(IQAir, 16–17 Apr)",
            ha="center", fontsize=8)
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b"))
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8, frameon=False)
    ax.set_title("Case study April 2020 (Chornobyl Exclusion Zone / Zhytomyr fires): upwind fire load peaks,\n"
                 "but CAMS PM2.5 stays low — CAMS misses the smoke, so it cannot validate the fire→PM link")
    save(fig, "07_kyiv_case_april2020.png")


def kyiv_sources():
    t = pd.read_csv(os.path.join(PROC, "kyiv", "kyiv_highpm_load_by_type.csv"), index_col=0).iloc[:, 0]
    r = pd.read_csv(os.path.join(PROC, "kyiv", "kyiv_highpm_load_by_region.csv"), index_col=0).iloc[:, 0].head(8)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    t = t.reindex([x for x in lib.TYPE_ORDER if x in t.index])
    axes[0].barh(t.index[::-1], t.values[::-1] * 100, color=[lib.TYPE_COLORS[x] for x in t.index[::-1]])
    axes[0].set_xlabel("% of upwind fire load")
    axes[0].set_title("Burn types")
    axes[1].barh(r.index[::-1], r.values[::-1] * 100, color="#21918c")
    axes[1].set_xlabel("% of upwind fire load")
    axes[1].set_title("Source regions")
    fig.suptitle("What burns upwind of Kyiv on poor-air days with high upwind load (2013–2026; association, "
                 "not causation)", fontsize=11)
    save(fig, "08_kyiv_upwind_sources.png")


if __name__ == "__main__":
    validation_table()
    w = ukraine_weekly()
    calendar_heatmap(w)
    critical_weeks(w)
    burn_type_calendar()
    kyiv_case()
    kyiv_sources()
