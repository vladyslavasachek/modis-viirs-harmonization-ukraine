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

## Pre-Stage-5 fixes (19:08)
- **Trend area check:** the open-ended conflict oblasts were already excluded in all years (100% flagged before and after 2022). But the north Feb–Apr 2022 window removed 4,733 of 7,344 VIIRS-equivalent hot spots in the 4 northern oblasts in 2022 alone. Trends now use a strict **constant area**: open-ended oblasts excluded in all years, north window kept. Computed from detection weights.
- **Corrected trend table** (`trend_robustness_ukraine.csv`):

  | Domain | Series | Period | Trend | p |
  |---|---|---|---|---|
  | Constant area | Harmonized | 2001–25 | −1,914/yr (−4.2%/yr) | 0.009 |
  | Constant area | Raw Aqua | 2003–25 | −3.5%/yr | 0.064 |
  | Constant area | Raw S-NPP | 2012–25 | −11.2%/yr | 0.019 |
  | Constant area | Harmonized | 2001–21 | −2.4%/yr | **0.32** |
  | All Ukraine | Harmonized | 2001–21 | −2.8%/yr | 0.20 |
  | All Ukraine | Raw Aqua | 2003–21 | −0.5%/yr | 0.84 |
  | All Ukraine | Raw S-NPP | 2012–21 | −4.0%/yr | 0.60 |
  | All Ukraine incl. conflict | Harmonized | 2001–25 | −2.8%/yr | 0.08 |

  - The signs agree everywhere. **The significant constant-area decline depends on 2022–2025**; there is no significant pre-invasion trend.
- **Baseline toggle** (sidebar; applies to KPI, climatology, anomalies, critical periods, season split): 2001–2025 vs 2012–2025 VIIRS-only. The default is VIIRS-only for AOIs < 20 cells of 0.1°.
- Limitations now cover: point estimates in percentiles, weekly resolution, zero-week ranks, fixed season split, the pooled threshold, and the trend's dependence on 2022–25.
- Belarus FIRMS arrived:
  - S-NPP: 86,249 archive (2012-01-26 → 2026-07-31) + 940 NRT;
  - NOAA-20: 41,311 archive (2018-04-03 → 2026-07-31) + 970 NRT.

## Stage 5 — Burn-type classification, full version (done)
- `src/events.py`: space-time clustering of VIIRS **S-NPP** detections (2012 → 2 Nov 2026), then NOAA-20 NRT after S-NPP ends.
  - Link rule: ≤ 1 km and ≤ 2 days (ST-DBSCAN, min_samples = 1 = connected components of the space-time neighbour graph); clustered per calendar year.
  - Justification: 1 km links adjacent 375 m pixels incl. diagonals and off-nadir growth; 2 days bridges one missed/clouded overpass. S-NPP only, so event sizes stay comparable over time.
  - Result: 1.03 M detections → **446,676 events** (193,553 with ≥ 2 detections). Median 1 detection, 99th percentile 14, max 1,612. Runs in 4 s.
- Features per event: duration, n_det, mean/max/sum FRP, night fraction, extent, land-cover shares, day of year, recurrence (1 km cell burned in the previous 5 years, MODIS + VIIRS history), conflict fraction.
- **Rule types** (events / detections):

  | Type | Events | Detections |
  |---|---|---|
  | Crop residue | 258,596 | 589 k |
  | Grass | 95,626 | 201 k |
  | Forest | 66,737 | 148 k |
  | Other | 25,502 | 84 k |
  | Peat / wetland | 215 | 13 k |

  - Peat events: median 20 detections, 4.5 days, 74% at night, about 3 MW.
  - Share of harmonized activity (all years, both domains): crop 65%, grass 16%, forest 12.5%, other 4.4%, peat 1.5%.
- **GMM cross-check** (k = 5, standardized features): agrees with the rules for **85%** of events; ARI 0.45, NMI 0.51. Four clusters align with forest / spring crop / autumn crop / grass; the fifth = large long mixed events.
- **Known events** (`events_known_cases.csv`):
  - CEZ April 2020: **98.5%** of detections typed forest (rest grass/other). Consistent with the reported forest/grass fires; 365 km² burned per Newman-Thacker & Turnbull 2021.
  - Kyiv-region peat fires, Aug–Oct 2015 (Vyshhorod; Babel, EcoPolitic, Karazin University journal): only **8%** typed peat, 63% grass. Drained peatlands appear as grassland in WorldCover, and the VIIRS signal was short, daytime and intense. **The peat rule has low recall for drained peat.**
- `validation_events_worldview.csv`: 40 events, 8 per type (n_det ≥ 3, 2012–2025, non-conflict). Worldview links "during" and "after (burn scar)" use a MODIS true-colour base plus S-NPP, NOAA-20 and MODIS fire overlays. S-NPP imagery layers are flagged in Worldview after the end of S-NPP. Columns `my_label`, `my_notes` are empty.
- Per-detection `burn_type` added to `detections.parquet`: S-NPP from its event (100% matched); MODIS / NOAA-20 / pre-2012 from land cover only (`type_from_event = False`, lower confidence).
- App page "Burn types":
  - monthly calendar by type;
  - shares by year and season;
  - event map (Okabe-Ito colours);
  - "How does this fire compare?" (percentiles of size/extent/total FRP/peak FRP/duration among same-oblast events);
  - validation expander.
- `data/processed/` = 39 MB.

## Stage 6 — Kyiv air (done, 19:45)
- `src/extra_regions.py` reads any `data/raw/firms_*/` folder with no code change.
  - Belarus: 129,470 raw → 89,205 filtered detections; 26,077 events.
  - Detections by type: grass 36%, forest 24%, crop 23%, other 12%, peat 4%.
- `src/kyiv_air.py`, daily upwind fire load (UFL):
  - Inputs: S-NPP stream within 300 km of Kyiv; bearing within ±30° of the ERA5 100 m wind-from direction; FRP × exp(−d/100 km).
  - Lags 0/1/2 with weights 1/0.5/0.25, using the mean transport-window wind.
  - Split by type and region; conflict-flagged detections excluded.
  - 340 k detections within 300 km; Belarus ≈ 3.6% of distance-weighted FRP.
- PM: OpenAQ needs a key (401 without), so the **CAMS** series via Open-Meteo is used: complete hourly 2013-01-01 → 2026-11-14. Circularity caveat stated.
- **Regression** (log PM2.5 minus trailing level; controls: month, wind speed/dir, BLH, precip, Tmax, RH):
  - Chronological test 2022–26: meteorology-only R² 0.651 vs + UFL 0.650. **No added skill.**
  - Full-sample HAC: lag-1 UFL +0.67% PM per doubling of (1+UFL), p = 0.009; lag 0/2 not significant.
- **High-PM days** (> p90 = 23.1 µg/m³): 12.2% have high UFL vs 6.7% of all days (lift 1.8). Spring 30% vs 12%; Jul–Oct 36% vs 9%.
  - Load on those days by type: crop 65%, grass 15%, forest 10.5%, other 9%, peat 0.2%.
  - By region: Kyiv obl 34%, Cherkasy 15%, Zhytomyr 12%, Kirovohrad 10%, Vinnytsia 9%.
- **Key diagnostic, April 2020:** UFL among the highest in the record (10 Apr 4,129; 17 Apr 3,248; NW wind from CEZ/Zhytomyr). Kyiv was the world's most polluted city on 16–17 Apr 2020 (IQAir AQI 238–380; Hromadske, Böll Stiftung), but **CAMS PM2.5 at Kyiv = 6–16 µg/m³**. CAMS misses the episode, so CAMS cannot validate the fire–air link. **Ground-station PM is needed (OpenAQ key / Kyiv city network).**
- Case studies in the app: April 2020, Aug–Sep 2015, Oct 2019, Mar 2015.
- App page "Kyiv air": UFL calendar, by-type and by-region charts, high-PM statistics, case studies, regression with caveats, CSV download.

## Peat improvement with PEATMAP (time-boxed, done 19:50)
- PEATMAP Europe (Xu et al. 2018, doi:10.5518/252, CC BY 4.0) downloaded, 170 MB, raw only; the `Other_European_Peatland` layer was used. `src/peat.py`.
- New rule: an event is peat if ≥ 50% of its detections are on mapped peat (regardless of WorldCover), OR it meets the previous wetland + smouldering signature.
- Effect:
  - peat events 215 → **6,856**; detections 13.4 k → 34.7 k (from crop 3.6 k, forest 4.1 k, grass 5.9 k, other 7.7 k);
  - share of harmonized activity typed peat: 1.5% → 2.5%;
  - Belarus peat share 4% → 16% (Polesia).
- **Before/after, Kyiv-region 2015 check: unchanged, 7.8% peat.** PEATMAP has only 8 polygons in 29–32E / 50–51.6N, and the nearest is ~70 km from the Vyshhorod fires (83 km from Kyiv). The limitation stays for the Kyiv region.
- CEZ April 2020 unchanged (98.5% forest).
- Kyiv high-PM/high-load days: peat share of load 0.2% → 3.2%.
- `validation_events_worldview.csv` is protected from overwriting; event IDs are unchanged.
