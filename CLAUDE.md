# Project: Harmonized burning-activity calendar for Ukraine, with a focus on Kyiv's air
NASA Space Apps Challenge 2026 — "Harmonization of MODIS and VIIRS Hot Spots"

## 0. How to work (read first)
- Time is very limited (submission is due tonight). Work strictly in the stage order below.
  A simple version that works and is honest beats a complex one that is unfinished.
- After each stage: (1) leave the code in a working state, (2) make a git commit with a clear
  message, (3) append a short entry to PROGRESS.md (what was done, key numbers, open issues),
  and (4) give me a short report with the key numbers and any figure paths.
- If a stage is taking too long, simplify it (use the fallback named in the stage), tell me
  why, and move on.
- Before writing the ETL, inspect the raw files (names, columns, row counts, date ranges per
  sensor) and report what you found.
- Use real data only. Never invent numbers or substitute synthetic data for real results.
  If a download fails or a source does not have the coverage you need, stop and tell me.
- Python 3.10+, pandas/pyarrow, numpy, scipy, statsmodels, scikit-learn, geopandas/shapely,
  rasterio, plotly, Streamlit. Create a virtual environment and requirements.txt.
- Keep secrets in .env (never commit it). Raw data goes in data/raw/ (gitignored).
  Processed tables that the app needs go in data/processed/. Keep them small (target <50 MB
  total, parquet) so they can be committed and the app can be deployed.

## 1. Context
Satellites have recorded active fire hot spots for over 20 years, but the record is split
across sensors whose data cannot be compared directly:
- MODIS (Terra from 1 Nov 2000, Aqua from 4 Jul 2002): ~1 km pixels. The FRP detection limit
  is several MW, and higher off-nadir and in daytime.
- VIIRS 375 m (S-NPP 20 Jan 2012 – 2 Nov 2026; NOAA-20 from 1 Apr 2018; NOAA-21 from 2023):
  detects small and weak fires down to ~0.5 MW, and has much better night-time performance.
  One fire that MODIS sees as one pixel is often several VIIRS pixels.
  As a result, a naive merged record jumps several-fold in 2012. That is a sensor artefact,
  not more fires.
- S-NPP data delivery has ended, so NOAA-20 must carry the reference record forward.

Official goal: a web application that harmonizes these records into a burning-activity
calendar based on active fire hot spots. It must let emergency responders, scientists and
land managers examine historical fire patterns, unusual conditions and critical periods
within a selected area of interest (AOI).

Our angle: Kyiv has recurring seasonal smog episodes from uncontrolled burning (peat,
crop-residue, grass and forest fires). We fully solve the official task for all of Ukraine,
and add two things that other teams do not have:
(A) automatic classification of burning type, and
(B) a link between burning and Kyiv's air quality.

## 2. Scientific rules (must never be violated)
1. A hot spot detection is an observation, not a fire. Never add MODIS and VIIRS counts
   together as "number of fires".
2. No detection does not mean no fire (clouds, overpass gaps, detection limits).
3. Harmonized values are model estimates with uncertainty. Always show intervals.
4. Every model is compared against a simple baseline using chronological / leave-one-year-out
   validation.
5. Descriptive and statistical analysis only. No ignition forecasts.
6. Associations are not causation. Say so wherever it matters.
7. Burning-type labels are model output, not ground truth. Show confidence and say how they
   were validated.

## 3. Data
Already downloaded (inspect first): data/raw/firms/. These are FIRMS archive downloads for
Ukraine (CSV): MODIS C6.1, VIIRS S-NPP 375 m, VIIRS NOAA-20 375 m, possibly NOAA-21.
FIRMS deliveries usually contain both `fire_archive_*` (standard processing) and `fire_nrt_*`
(near-real-time, recent months) files.
- Keep a `processing` column (standard / nrt) and do not mix the two silently.
- NRT files may lack the `type` column. Handle that, for example by applying a
  static-source mask built from the standard data.
- Filters:
  - `type == 0` (presumed vegetation fire). Exclude type 2 (static land sources such as
    industry) and build a persistent-hot-spot mask from them.
  - MODIS confidence ≥ 30.
  - VIIRS confidence in {n, h}.
  - Run a short sensitivity check on these thresholds.

Download yourself (cache everything in data/raw/):
- Boundaries: geoBoundaries Ukraine ADM0 + ADM1. Chornobyl Exclusion Zone: use a published
  polygon if you can find one, otherwise a documented approximate polygon.
- Land cover: ESA WorldCover 10 m 2021 (public AWS bucket). Classes used: cropland, tree
  cover, grassland, shrubland, herbaceous wetland, built-up, bare.
  Fallback if too heavy: sample WorldCover only at detection locations.
- Weather: Open-Meteo Historical Weather API (ERA5, no key). Daily and hourly as needed:
  temperature_2m_max, precipitation_sum, et0_fao_evapotranspiration, wind_speed_10m,
  wind_direction_10m (and 100 m wind if available), relative humidity.
- Air quality for Kyiv, in this order of preference:
  1. ground-station PM2.5/PM10 (OpenAQ; free API key, ask me for it if needed),
  2. Open-Meteo Air Quality API (CAMS).
  Report the actual historical coverage before using it.
  IMPORTANT: CAMS assimilates fire emissions derived from MODIS/VIIRS FRP (GFAS), so a
  correlation between CAMS PM and FIRMS fires is partly circular. Prefer ground stations,
  and state this caveat explicitly if you use CAMS.
- Optional: Meteomatics API (credentials in .env as METEOMATICS_USER / METEOMATICS_PASSWORD).
  Use it only for parameters that Open-Meteo lacks. Check its documentation and quota first.

Attribution (in README and app): NASA FIRMS; MODIS C6.1 active fire (MCD14ML / MCD14DL);
VIIRS 375 m active fire (VNP14IMG, VJ114IMG); ESA WorldCover; ERA5 via Open-Meteo; any
air-quality source used.

## 4. Stages

### Stage 1 — ETL (required)
- One clean hot spot table (parquet) with columns: sensor, satellite, datetime UTC, day/night,
  lat, lon, FRP, scan, track, confidence, type, processing.
- Base aggregate: 0.1° cell × UTC day × day/night × sensor (count, FRP sum, FRP max).
- Attach WorldCover class to every detection (and the class shares to every 0.1° cell).
- Report row counts and date ranges per sensor, plus a raw annual-count plot per sensor.

### Stage 2 — Harmonization (required, core of the challenge)
Reference scale: VIIRS S-NPP. Goal: estimate VIIRS-equivalent activity from MODIS for
2000–2012 (and for S-NPP outage days). Then continue the record with NOAA-20 after S-NPP ends.
- Train on the overlap period, 2012–2025.
- Compare these models, working at the cell × week level (0.25° or 0.5°, whichever validates
  better):
  - B0: a single global VIIRS/MODIS ratio.
  - B1: ratio by month × day/night.
  - M1: Negative Binomial GLM of the VIIRS count on log1p(MODIS count), log1p(MODIS FRP),
    month, day/night and land-cover shares, with an explicit term for cells where MODIS = 0
    (the "missed" small fires only VIIRS sees). Use partial pooling / empirical-Bayes
    shrinkage for small strata.
  - Fallback if M1 is slow: B1 plus a land-cover-specific ratio.
- Terra-only period (Nov 2000 – Jul 2002): either a separate Terra-only sub-model or an
  explicit "lower confidence" flag.
- Validation: leave-one-year-out on 2012–2025, predicting VIIRS from MODIS only.
  - Metrics: WAPE and bias at Ukraine × month, Ukraine × year, oblast × month and
    cell × month.
  - Output: a results table B0 vs B1 vs M1. Choose the model by validation, not complexity.
- Uncertainty: bootstrap over training years, giving 80% and 95% intervals.
- Continuity: calibrate NOAA-20 to S-NPP over 2018–2026 (ratio by month × day/night, with
  validation), then use calibrated NOAA-20 after S-NPP ends.
- Key figure: raw MODIS, raw VIIRS and the harmonized series on one axis for Ukraine. The
  2012 jump must disappear.

### Stage 3 — Calendar, anomalies, critical periods (required)
For any AOI (computed from the cell table):
- Heatmaps of week-of-year × year and month × year of harmonized activity, with uncertainty
  available on hover.
- Climatology per week: median and 10–90% range, computed from baseline years that exclude
  the year being evaluated.
- Anomaly: the percentile of a week or month relative to the baseline. Do not use z-scores;
  the data are skewed and zero-inflated.
- Critical periods: for each week of the year, the probability that activity exceeds the
  baseline 75th percentile (bootstrap interval). Highlight weeks above 50%.
- Season metrics per year: onset (10% of annual cumulative), peak and end (90%), computed
  separately for spring and autumn if the season is bimodal.

### Stage 4 — Minimal web app (required; build early, extend later)
- Streamlit app (app/app.py) with an English UI.
- AOI selection:
  - presets: Ukraine, "Kyiv airshed" (Kyiv city + Kyiv, Zhytomyr and Chernihiv oblasts),
    the Chornobyl Exclusion Zone, and every oblast;
  - a drawn polygon (streamlit-folium Draw).
- Date range selection.
- Pages at this stage: Calendar; Anomalies & critical periods; Harmonization (raw vs
  harmonized, validation table); Data, methods & limitations.
- CSV download of the calendar for the selected AOI.
- The app must run fast using only data/processed/.

### Stage 5 — Original part A: burning-type classification
Unit: fire events. Build them with a space–time clustering of VIIRS detections (2012+), e.g.
ST-DBSCAN with ~1 km / 1–2 day linkage; document and justify the parameters.
- Features per event: dominant WorldCover class and class shares; duration (days); number of
  detections; mean/max/total FRP; fraction of night detections; spatial extent; day of year;
  recurrence of the same cells in earlier years; distance to settlements (optional).
- Classes (transparent rules first):
  - crop-residue / field burning: cropland dominant, short, spring/late-summer/autumn peaks;
  - grass / abandoned-land burning: grass or shrubland dominant;
  - forest fire: tree cover dominant;
  - peat / wetland smouldering: herbaceous wetland or known peat areas, long duration, high
    night fraction, low FRP per detection;
  - other / uncertain.
- Then cross-check with unsupervised clustering (GMM or k-means on standardized features).
  Report how well the clusters agree with the rule labels.
- Validation:
  1. known events, including the April 2020 Chornobyl Exclusion Zone fires (forest/grass) and
     any well-documented peat fire episodes near Kyiv you can find in reliable sources
     (cite them);
  2. a CSV of 40 random events with NASA Worldview links (date + coordinates), so that I can
     check them by eye. Leave a column for my label.
- Pre-2012 (MODIS only): assign land-cover-based type shares only, flagged as lower
  confidence.
- App page "Burn types": the calendar split by type for the AOI, a map of events by type, and
  type shares by season and year.
- Fallback if time is short: land-cover-only typing at detection level, plus the night
  fraction for the peat flag.

### Stage 6 — Original part B: burning and Kyiv's air
Question: which burning, and from where, coincides with poor air in Kyiv?
- Daily "upwind fire load" for Kyiv:
  - sum the FRP of detections within 300 km whose bearing from Kyiv lies within ±30° of the
    daily upwind direction;
  - weight by distance decay;
  - use the previous 1–2 days as well, to allow for transport;
  - split the load by burning type and by oblast.
  Use ERA5 wind (100 m if available, else 10 m).
- Relate it to Kyiv PM2.5 (or PM10). Use a regression of the daily PM anomaly on upwind fire
  load, controlling for season, wind speed, precipitation and temperature. Compare it to a
  meteorology-only baseline on a chronological split.
- Outputs:
  - the share of high-PM days (above the 90th percentile) with high upwind fire load;
  - which burning types and oblasts dominate on those days;
  - 2–3 case studies (e.g. April 2020, if the data cover it).
  Frame as "where and what to control to protect Kyiv's air", stated as an association.
- App page "Kyiv air".
- Fallback if no usable PM data: show only the upwind fire load calendar for Kyiv by type and
  source region, and state that air-quality validation is future work.

### Stage 7 — Optional, only if time remains
- Season-shift trends: Mann–Kendall + Sen's slope on onset/peak/end dates, 2003–2025. Check
  robustness on raw Aqua-MODIS-only and on VIIRS-only series; timing is relative, so
  agreement across all three means the result is not a harmonization artefact.
- Weather drivers: regression of monthly harmonized activity anomalies on weather anomalies
  (temperature, precipitation, cumulative P − ET0 over 30/60/90 days). Residuals = the
  "non-weather" anomaly.

### Stage 8 — Deliverables (reserve time for this)
- README.md (English): problem, data, method, validation numbers, burn-type and Kyiv-air
  results, limitations, how to run, attribution.
- figures/: 6–8 PNGs for slides (the 2012 jump removed, validation table, calendar,
  critical weeks, burn-type calendar, Kyiv upwind load vs PM, a case study).
- SUBMISSION.md: a project summary of ≤ 200 words; a "what it does / how it works /
  benefits / limitations / tools" section; the list of data sources; and a note on the use of
  AI tools.
- Make sure the app runs from a clean clone using only committed processed data, and is
  ready for deployment to Streamlit Community Cloud.
