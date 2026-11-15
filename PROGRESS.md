# PROGRESS

## Setup (2026-11-15 17:45 Kyiv)
- git repo, .venv (Python 3.13), requirements.txt, .gitignore (.env, data/raw/ excluded).

## Stage 1 — ETL (done)
- Scripts: `src/etl.py` → `src/worldcover.py` → `src/aggregate.py` (~75 s total).
- After filters (type 0 or NRT static mask; MODIS conf ≥30; VIIRS n/h; inside UKR ADM0), detections are:
  - MODIS 626,169 standard + 2,874 NRT (2000-11-01 → 2026-11-13);
  - S-NPP 994,530 + 40,299 (2012-01-25 → 2026-11-01);
  - NOAA-20 468,087 + 42,759 (2018-04-01 → 2026-11-13);
  - NOAA-21 182,525 + 37,760 (2023-03-01 → 2026-11-14).
- NRT static mask: 1,540 cells (0.004°, ≥3 type-2 hits). On standard data it captures 96.8% of type 2 and wrongly masks 0.95% of type 0.
- The FIRMS download is clipped to Ukraine, with almost no detections outside it. **The Stage 6 cross-border load (Belarus/Polesia) needs an extra FIRMS download.**
- Conflict flag (oblast × date window): 590,208 detections in Ukraine flagged.
- WorldCover 2021 v200 is read from the COG overview ×8 (~67 m) and attached to each detection, with 0.1° cell shares. Land-cover mix:
  - MODIS: crop 68%, grass 16%, tree 13%, wetland 2%;
  - S-NPP: crop 56%, grass 23%, tree 15%, wetland 4%.
  - The overview ×8 smooths small patches.
- MODIS 6.03 → 61.03 check: the S-NPP/Aqua ratio is 5.74 for 2017–22 vs 5.45 for 2023–25 (inside the 2017–22 range). There is no count step, but Aqua median FRP rose from about 16–17 to about 21 MW. Use count-based harmonization and treat FRP features with care after 2023.
- Outputs:
  - `data/processed/agg_cell_day.parquet` (8.5 MB, 1.25 M rows);
  - `cell01_landcover_shares.parquet`;
  - `qc_*.csv`;
  - `figures/01_raw_annual_counts.png`.
- Open issues: there are no Belarus detections; WorldCover is approximate (overview).

## Stage 2 — Harmonization (done, 18:27)
- Scripts: `src/harmonize.py` (library), `stage2_validate.py` (LOYO, B0/B1/M1, 0.25° and 0.5°), `stage2_validate_extra.py` (B2 = B1 × land cover, M1 without FRP), `stage2_noaa20.py`, `stage2_final.py`.
- Setup: cell × week × day/night; target S-NPP count; LOYO over 2012–2021 (pre-invasion); 65 S-NPP outage days (544 outage weeks incl. the pre-launch weeks) excluded.
- M1 = NB2 GLM by our own IRLS (identical to statsmodels to 1e-6, 15× faster); α = 0.21 (Terra+Aqua, 0.5°).
- **LOYO results at 0.5° (Terra+Aqua):**

  | Model | Ukraine × month WAPE | Ukraine × year WAPE | Oblast × month WAPE | Cell × month WAPE |
  |---|---|---|---|---|
  | B0 | 16.2% | **4.8%** | 21.7% | 37.0% |
  | **B1** | **9.7%** | 6.4% | **17.1%** | **34.7%** |
  | M1 | 14.7% | 10.0% | 23.8% | 37.5% |
  | B2 | 10.1% | 5.9% | 17.3% | 35.3% |

  - The 0.25° grid is no better for any model (M1 is worse).
  - **Chosen: B1 at 0.5°** (shrunk month × day/night ratio, K = 200). The more complex M1 does not validate better.
  - Terra-only sub-model (pre Jul 2002): B1 gives 11.9% / 7.1% / 21.6% / 45.6%.
- **Intervals:** the year-bootstrap of ratios alone is far too narrow (80% interval covers 14% at Ukraine × month). Replaced by LOYO-residual-calibrated predictive intervals, binned by predicted size and by AOI scale (cell / oblast / Ukraine). Cross-validated coverage: 80% → 72–81%, 95% → 90–96% across all levels.
- **Sensitivity (fit 2012–2025, conflict detections excluded):** 22/24 ratios fall inside the main 95% bootstrap range. The two outside are May and Sep night (−20%). The main fit applied out-of-sample to 2022–25 gives bias +5.8%, Ukraine × year WAPE 8.0%.
- **NOAA-20 → S-NPP:** a single ratio of 1.011 (95% 1.000–1.023). LOYO Ukraine × month WAPE is 3.7%; month × day/night ratios do not help.
- **Conflict "excl" variant:** the open-ended conflict oblasts are now removed in **all** years (domain-consistent); north Feb–Apr 2022 only in that window. Otherwise the series shows an artificial 2022 drop, since those oblasts held 65% of 2022 detections.
- Outputs:
  - `harmonized_cells.parquet` (1.2 MB);
  - `harmonization_params.npz`;
  - `interval_calibration_*.csv`;
  - `ukraine_*_harmonized*`;
  - `validation_*.csv`;
  - `figures/02_harmonized_vs_raw_ukraine.png`.
- Open issues:
  - B1 is zero wherever MODIS saw nothing, so the totals are calibrated but the spatial pattern of VIIRS-only small fires before 2012 is missing (cell × month WAPE 35%).
  - Aqua FRP drift after 2023 does not affect B1, which uses counts only.

## Stage 4 — Minimal app (done)
- `src/build_app_data.py` builds `data/processed/app/`:
  - `weekly_cells.parquet` (0.1° × week, 2.6 MB): B1 is a cell-independent ratio, so it is applied at 0.1°. It reproduces the Stage 2 annual totals to 1e-7.
  - `cells.parquet` (AOI membership).
  - `boundaries.geojson`.
  - `detections.parquet` (2.40 M detections, 29 MB) with a per-detection **harmonization weight**. The sum of weights over any area or time window equals the harmonized value; this drives the daily month-to-date KPI and the exact circle/polygon AOIs.
- CEZ polygon: OSM relation 3311547 (2,627 km², official ≈ 2,600 km²).
- `app/app.py` + `app/lib.py`. Pages: Calendar; Anomalies (placeholder until Stage 3); Harmonization; Data, methods & limitations.
- AOI: presets (Ukraine, Kyiv airshed, CEZ, 27 oblasts), point + radius, drawn polygon.
- Domain toggle, default "Ukraine without conflict-affected oblasts, all years".
- Calendar page:
  - KPI "this month to date vs 2001–2025 median" with intervals and percentile; for Ukraine, Nov 1–14 2026 = 44 vs median 264, the 20th percentile;
  - weekly series with 80/95% bands;
  - week × year and month × year heatmaps (viridis, intervals on hover);
  - annual bars with 95% error bars and Sen's slope + Mann–Kendall (Ukraine excl. domain 2001–2025: −1,918/yr, 95% CI −3,336 to −622, p = 0.009);
  - top 10 months and years;
  - calendar CSV and hot-spot CSV downloads; PNG via the plot camera icon.
- Tested in the browser: preset, point + radius and drawn rectangle. No server errors.
- Total `data/processed/` is 44 MB.

## Checks before Stage 3 (18:50)
- **Trend robustness** (`data/processed/trend_robustness_ukraine.csv`). Sen's slope + Mann–Kendall, Ukraine without conflict-affected oblasts:

  | Series | Period | Trend | p |
  |---|---|---|---|
  | Harmonized | 2001–25 | −1,918/yr (−4.3%/yr) | 0.009 |
  | Raw Aqua-only | 2003–25 | −280/yr (−3.5%/yr) | 0.07 |
  | Raw S-NPP-only | 2012–25 | −3,611/yr (−11.3%/yr) | 0.02 |

  - **All three agree in sign** (decline).
  - All of Ukraine (incl. conflict oblasts): all negative, none significant (p 0.08–0.35).
  - Shown on the Calendar page for any AOI. README must list all three.
- **Month-to-date KPI:** confirmed the same 1–14 Nov window in each baseline year (independent recomputation, max diff 0). The panel states the current month is NRT (NOAA-20 calibrated).
- **Size:** `data/processed/` went from 44 MB to **33 MB**. Detections are quantized (lat/lon 1e-4°, FRP 0.1 MW, minutes) from 29.5 to 19.9 MB; zstd level 19; compact aggregate dtypes.

## Stage 3 — Calendar, anomalies, critical periods (done, 19:00)
- `src/calendar_stats.py`: leave-one-year-out climatology (median, 10–90%), mid-rank percentile anomalies (weekly/monthly), critical periods, bimodal season split and onset/peak/end (10/50/90% cumulative).
- **Deviation from the spec (documented in the app):** the critical-period threshold is the 75th percentile *pooled over all weeks* of the baseline years (without the evaluated year). A per-week threshold would be exceeded in ~25% of years by construction. 95% CI by bootstrap over years.
- Ukraine (excl. domain):
  - critical weeks 12–15 (late Mar–mid Apr) and 29–34 + 38 (late Jul–Sep);
  - bimodal, split at week 23;
  - median spring onset/peak/end at day 74 / 88 / 116;
  - median summer–autumn onset/peak/end at day 200 / 228 / 277.
- Kyiv airshed: April 2020 = **100th percentile**; critical weeks 12–16 and 31–38.
- App page "Anomalies & critical periods":
  - climatology vs selected year;
  - RdBu percentile heatmap (centred at 50);
  - critical-weeks bars with CI;
  - season timing;
  - CSV downloads.
