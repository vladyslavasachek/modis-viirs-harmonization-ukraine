"""Stage 2 validation: LOYO 2012-2021, B0/B1/M1, grids 0.25 and 0.5, Terra+Aqua and Terra-only;
sensitivity fit 2012-2025 with conflict-flagged detections excluded."""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from harmonize import (B0, B1, M1, LC, PROC, SENS_YEARS, TRAIN_YEARS, build_table, cell_id, full_grid,
                       load_det, loyo, lc_shares, metrics, outage_weeks, prepare)


def cell_oblast_map(d, res):
    c = d.assign(cell=cell_id(d["lat"], d["lon"], res)).groupby(["cell", "oblast"]).size().reset_index(name="k")
    return c.sort_values("k").drop_duplicates("cell", keep="last")[["cell", "oblast"]]


def run_variant(d, res, sats, years, outage, exclude_conflict=False, models=("B0", "B1", "M1")):
    t = time.time()
    cells = np.unique(cell_id(d["lat"], d["lon"], res))
    m, v = build_table(d, res, years, modis_sats=sats, exclude_conflict=exclude_conflict)
    g = prepare(full_grid(cells, years, m, v), lc_shares(res), cell_oblast_map(d, res), outage)
    g = g[~g["outage"]].copy()
    alpha = M1().fit(g).alpha
    fac = {"B0": B0, "B1": B1, "M1": lambda: M1(alpha=alpha)}
    preds = loyo(g, years, {k: fac[k] for k in models})
    rows = []
    for k, p in preds.items():
        r = metrics(p)
        r.update(model=k, grid=res, modis="+".join(sats), years=f"{years[0]}-{years[-1]}",
                 conflict_excluded=exclude_conflict)
        rows.append(r)
    print(f"variant res={res} sats={sats} n_rows={len(g):,} cells={len(cells)} alpha={alpha:.3f} "
          f"time={time.time() - t:.0f}s", flush=True)
    return pd.DataFrame(rows), g, preds, alpha


def main():
    d = load_det()
    outage, od = outage_weeks(d)
    print("S-NPP outage days (>=2012-01-25):", int((od >= pd.Timestamp('2012-01-25')).sum()),
          "outage weeks:", len(outage), flush=True)
    res_all = []
    for res in (0.5, 0.25):
        for sats in (("Terra", "Aqua"), ("Terra",)):
            r, g, preds, alpha = run_variant(d, res, sats, TRAIN_YEARS, outage)
            res_all.append(r)
            if res == 0.5 and sats == ("Terra", "Aqua"):
                # save LOYO predictions at Ukraine x week for figures
                lp = pd.concat([p.groupby(["year", "woy"])[["pred", "v_n", "m_n"]].sum().assign(model=k)
                                for k, p in preds.items()]).reset_index()
                lp.to_parquet(os.path.join(PROC, "loyo_ukr_week.parquet"), index=False)
    tab = pd.concat(res_all, ignore_index=True)
    tab.to_csv(os.path.join(PROC, "validation_loyo.csv"), index=False)
    cols = ["modis", "grid", "model", "ukr_month_WAPE", "ukr_month_bias", "ukr_year_WAPE", "ukr_year_bias",
            "oblast_month_WAPE", "cell_month_WAPE"]
    print(tab[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
