"""NOAA-20 -> S-NPP continuity calibration (ratio by month x day/night), LOYO over 2018-2025."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from harmonize import B0, B1, PROC, build_table, cell_id, full_grid, load_det, loyo, metrics, outage_weeks, prepare, lc_shares
from stage2_validate import cell_oblast_map

YEARS = list(range(2018, 2026))
RES = 0.5


def main():
    d = load_det()
    outage, _ = outage_weeks(d)
    cells = np.unique(cell_id(d["lat"], d["lon"], RES))
    # predictor "m" = NOAA-20 counts, target = S-NPP counts (standard processing)
    m, v = build_table(d, RES, YEARS, modis_sats=("NOAA20",))
    g = prepare(full_grid(cells, YEARS, m, v), lc_shares(RES), cell_oblast_map(d, RES), outage)
    g = g[~g["outage"]]
    g = g[~((g["year"] == 2018) & (g["woy"] < 13))]  # NOAA-20 starts 2018-04-01
    preds = loyo(g, YEARS, {"N20_B0": B0, "N20_B1": B1})
    rows = []
    for k, p in preds.items():
        r = metrics(p)
        r["model"] = k
        rows.append(r)
    # identity (no calibration) as the trivial baseline
    ident = g.assign(pred=g["m_n"])
    r = metrics(ident)
    r["model"] = "N20_identity"
    rows.append(r)
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(PROC, "validation_noaa20_calibration.csv"), index=False)
    print(tab[["model", "ukr_month_WAPE", "ukr_month_bias", "ukr_year_WAPE", "oblast_month_WAPE", "cell_month_WAPE"]]
          .round(3).to_string(index=False))

    # Chosen by LOYO: single global ratio (B0); month x dn ratios do not validate better.
    final = B0().fit(g)
    rng = np.random.default_rng(1)
    yr = g.groupby("year")[["v_n", "m_n"]].sum()
    boots = []
    for b in range(1000):
        ys = rng.choice(YEARS, size=len(YEARS), replace=True)
        s = yr.loc[ys].sum()
        boots.append(s["v_n"] / s["m_n"])
    q = np.quantile(boots, [.025, .1, .9, .975])
    out = pd.DataFrame([{"model": "N20_B0", "ratio": final.r, "q025": q[0], "q10": q[1], "q90": q[2], "q975": q[3]}])
    out.to_csv(os.path.join(PROC, "noaa20_to_snpp_ratio.csv"), index=False)
    print(out.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
