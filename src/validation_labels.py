"""Score the manual Worldview labels (validation_events_worldview.csv) against the rule-based burn types.

- Confident labels = not "unclear". Accuracy on events whose predicted type is a definite class
  (crop residue / grass / forest / peat); "other / uncertain" makes no type claim and is reported separately.
- 95% Wilson interval. Scored for the labels as sampled (WorldCover-only rules, before PEATMAP) and for the current
  labels (with PEATMAP). The sample is stratified (8 per predicted type), not representative of type prevalence.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from harmonize import INTERIM, PROC, ROOT

MAP = {"field": "crop residue / field", "crop": "crop residue / field", "grass": "grass / abandoned land",
       "forest": "forest", "peat": "peat / wetland smouldering"}
DEFINITE = ["crop residue / field", "grass / abandoned land", "forest", "peat / wetland smouldering"]


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def main():
    path = os.path.join(ROOT, "validation_events_worldview.csv")
    with open(path) as f:
        sep = ";" if f.readline().count(";") > 5 else ","
    v = pd.read_csv(path, sep=sep)
    v["my_label"] = v["my_label"].astype(str).str.strip().str.lower()
    e = pd.read_parquet(os.path.join(INTERIM, "events_all.parquet"), columns=["event_id", "burn_type"])
    v = v.merge(e.rename(columns={"burn_type": "current_label"}), on="event_id", how="left")
    v["truth"] = v["my_label"].map(MAP)
    v["confident"] = v["my_label"] != "unclear"
    out, cms = [], {}
    for col, name in (("rule_label", "as sampled (WorldCover-only rules)"), ("current_label", "current (with PEATMAP)")):
        c = v[v["confident"] & v[col].astype(str).isin(DEFINITE)]
        k, n = int((c[col].astype(str) == c["truth"]).sum()), len(c)
        lo, hi = wilson(k, n)
        out.append(dict(labels=name, n_confident_total=int(v["confident"].sum()), n_scored=n, correct=k,
                        accuracy=k / n if n else np.nan, wilson95_lo=lo, wilson95_hi=hi))
        cms[name] = pd.crosstab(v.loc[v["confident"], col].astype(str), v.loc[v["confident"], "my_label"],
                                rownames=["predicted"], colnames=["manual label"])
    res = pd.DataFrame(out)
    unclear = v.groupby("rule_label")["confident"].agg(n="size", unclear=lambda s: int((~s).sum()))
    unclear["unclear_rate"] = unclear["unclear"] / unclear["n"]
    per_type = (v[v["confident"]].assign(ok=lambda x: x["rule_label"] == x["truth"])
                .groupby("rule_label")["ok"].agg(confident="size", correct="sum"))
    res.to_csv(os.path.join(PROC, "validation_manual_accuracy.csv"), index=False)
    unclear.join(per_type).to_csv(os.path.join(PROC, "validation_manual_by_type.csv"))
    for name, cm in cms.items():
        cm.to_csv(os.path.join(PROC, f"validation_manual_confusion_{'sampled' if 'sampled' in name else 'current'}.csv"))
    pd.set_option("display.width", 200)
    print(res.round(3).to_string(index=False))
    print(unclear.join(per_type).to_string())
    for name, cm in cms.items():
        print("\n", name); print(cm.to_string())
    print("\nlabel changes sampled->current:", int((v["rule_label"] != v["current_label"]).sum()))
    print(v.loc[v["rule_label"] != v["current_label"], ["event_id", "rule_label", "current_label", "my_label"]].to_string())


if __name__ == "__main__":
    main()
