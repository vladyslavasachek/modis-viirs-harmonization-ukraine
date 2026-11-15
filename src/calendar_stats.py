"""Stage 3: climatology, percentile anomalies, critical periods and season metrics for any AOI series.

Inputs are tidy frames with columns year, woy (0..51) or month (1..12), and point (harmonized value).
Baselines always exclude the year being evaluated (leave-one-year-out). Percentiles, not z-scores
(the data are skewed and zero-inflated).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

BASE_YEARS = list(range(2001, 2026))   # complete years


def _matrix(df: pd.DataFrame, period: str, years) -> pd.DataFrame:
    n = 52 if period == "woy" else 12
    cols = range(n) if period == "woy" else range(1, 13)
    return df.pivot_table(index="year", columns=period, values="point", aggfunc="sum").reindex(
        index=years, columns=cols).fillna(0.0)


def climatology(df: pd.DataFrame, period: str = "woy", exclude_year: int | None = None, years=BASE_YEARS):
    """Median and 10-90% range per week (or month) over baseline years, excluding `exclude_year`."""
    M = _matrix(df, period, [y for y in years if y != exclude_year])
    return pd.DataFrame({"median": M.median(), "p10": M.quantile(.1), "p25": M.quantile(.25),
                         "p75": M.quantile(.75), "p90": M.quantile(.9), "n_years": len(M)})


def percentile_anomalies(df: pd.DataFrame, period: str = "woy", years=BASE_YEARS, eval_years=None):
    """Mid-rank percentile of each (year, period) value among the OTHER baseline years (same period)."""
    eval_years = eval_years or sorted(df["year"].unique())
    M = _matrix(df, period, sorted(set(years) | set(eval_years)))
    B = M.loc[years].values
    out = {}
    for y in eval_years:
        base = B[[i for i, yy in enumerate(years) if yy != y]]
        v = M.loc[y].values
        out[y] = ((base < v).sum(0) + 0.5 * (base == v).sum(0)) / base.shape[0] * 100
    P = pd.DataFrame(out, index=M.columns).T
    P.index.name = "year"
    return P, M.loc[eval_years]


def critical_periods(df: pd.DataFrame, years=BASE_YEARS, q: float = 0.75, B: int = 1000, seed: int = 0):
    """For each week of year: probability that weekly activity exceeds the baseline q-quantile of ALL weekly
    values (pooled over weeks). The threshold for year y is computed without year y. 95% CI by bootstrap
    over years. (A per-week threshold would give ~25% for every week by construction.)"""
    M = _matrix(df, "woy", years)
    X = M.values
    exceed = np.zeros_like(X, dtype=bool)
    for i in range(len(years)):
        thr = np.quantile(np.delete(X, i, axis=0), q)
        exceed[i] = X[i] > thr
    p = exceed.mean(0)
    rng = np.random.default_rng(seed)
    bs = np.array([exceed[rng.integers(0, len(years), len(years))].mean(0) for _ in range(B)])
    thr_all = float(np.quantile(X, q))
    return pd.DataFrame({"woy": M.columns, "prob": p, "lo95": np.quantile(bs, .025, axis=0),
                         "hi95": np.quantile(bs, .975, axis=0), "critical": p > 0.5}), thr_all


def season_split(df: pd.DataFrame, years=BASE_YEARS):
    """Detect bimodality from the median climatology (3-week smoothed). Returns split week or None."""
    med = climatology(df, "woy", years=years)["median"].values
    sm = np.convolve(np.r_[med[-1], med, med[0]], np.ones(3) / 3, mode="valid")
    a, b = sm[5:22], sm[22:45]           # spring window ~Feb-May, summer/autumn window ~Jun-Nov
    if a.max() <= 0 or b.max() <= 0:
        return None
    pa, pb = 5 + int(a.argmax()), 22 + int(b.argmax())
    trough = pa + int(sm[pa:pb + 1].argmin())
    lower = min(sm[pa], sm[pb])
    return trough if sm[trough] < 0.7 * lower and lower >= 0.25 * max(sm[pa], sm[pb]) else None


def season_metrics(df: pd.DataFrame, years=None, split: int | None = None) -> pd.DataFrame:
    """Onset (10% of cumulative), peak (max week), end (90%) per year (and per season if split given).
    Weeks are reported as day-of-year of the week centre (woy*7+4)."""
    years = years or sorted(df["year"].unique())
    M = _matrix(df, "woy", years)
    segs = {"spring": (0, split), "summer-autumn": (split, 52)} if split is not None else {"year": (0, 52)}
    rows = []
    for y in years:
        v = M.loc[y].values
        for name, (a, b) in segs.items():
            s = v[a:b]
            tot = s.sum()
            if tot <= 0:
                rows.append(dict(year=y, season=name, total=0.0, onset_doy=np.nan, peak_doy=np.nan, end_doy=np.nan))
                continue
            c = np.cumsum(s) / tot
            on, en, pk = a + int(np.argmax(c >= 0.1)), a + int(np.argmax(c >= 0.9)), a + int(s.argmax())
            rows.append(dict(year=y, season=name, total=float(tot), onset_doy=on * 7 + 4, peak_doy=pk * 7 + 4,
                             end_doy=en * 7 + 4))
    out = pd.DataFrame(rows)
    out["length_days"] = out["end_doy"] - out["onset_doy"] + 7
    return out


def season_shift_trends(series: dict, split: int | None, periods: dict) -> pd.DataFrame:
    """Mann-Kendall + Sen's slope on onset / peak / end day-of-year per season.
    series: {name: weekly frame (year, woy, point)}; periods: {name: (y0, y1)}. Slopes in days per decade."""
    from scipy.stats import kendalltau, theilslopes
    rows = []
    for name, df in series.items():
        y0, y1 = periods[name]
        sm = season_metrics(df[df["year"].between(y0, y1)], years=list(range(y0, y1 + 1)), split=split)
        for (sea), g in sm.groupby("season"):
            g = g[g["total"] > 0]
            for m in ("onset_doy", "peak_doy", "end_doy"):
                if g[m].notna().sum() < 8:
                    continue
                s = theilslopes(g[m], g["year"])
                tau, p = kendalltau(g["year"], g[m])
                rows.append(dict(series=name, period=f"{y0}–{y1}", season=sea, metric=m.replace("_doy", ""),
                                 slope_days_per_decade=s.slope * 10, ci_lo=s.low_slope * 10, ci_hi=s.high_slope * 10,
                                 tau=tau, p=p, median_doy=g[m].median()))
    return pd.DataFrame(rows)
