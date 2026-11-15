# Submission — Harmonized burning-activity calendar for Ukraine

## Project summary (≤ 200 words)
Satellites have mapped fires for 25 years, but MODIS and VIIRS cannot be compared directly: VIIRS sees 5–6× more hot
spots, so merged records jump in 2012. We built a web app that harmonizes MODIS (2000–) to the VIIRS S-NPP scale,
continues the record with NOAA-20 after S-NPP ended, and turns it into a burning-activity calendar for any area of
Ukraine: preset, radius or drawn polygon. The model was chosen by leave-one-year-out validation (9.7% monthly error
for Ukraine), and every estimate carries a calibrated uncertainty interval. The app shows anomalies as percentiles,
critical weeks, season timing and long-term trends. It also downloads hot spots and calendars. Two additions go
beyond the challenge. Fire events are classified into crop-residue, grass, forest and peat burning (78% agreement in
a manual satellite-image check). Upwind burning is traced to Kyiv using wind data, including fires in Belarus. On
poor-air days, the upwind load comes mostly from crop-residue burning in Kyiv, Cherkasy and Zhytomyr oblasts. We
state openly what is not validated: model-based PM data miss the April 2020 smoke episode, and ground stations
exist only since 2023.

## What it does
- One consistent, uncertainty-aware fire-activity record for Ukraine, 2000–2026, at 0.1° × week, for any AOI and date range.
- Calendar heatmaps, weekly series with 80/95% intervals, "this month vs historical median" card, top months/years,
  trends; percentile anomalies, critical weeks, season onset/peak/end and season-shift trends.
- Burn-type calendar and event map; "How does this fire compare?" for any event.
- Kyiv air: which burning, and from where, is upwind of Kyiv on poor-air days.
- CSV downloads (calendar, individual hot spots, anomalies, season metrics, Kyiv daily table), PNG for every plot.

## How it works
- NASA FIRMS archive + NRT (MODIS C6.1, VIIRS S-NPP / NOAA-20 / NOAA-21), vegetation-fire filtering, static-source mask.
- Harmonization model chosen by leave-one-year-out validation (B0 / B1 / B2 / NB-GLM); B1 (ratio by month ×
  day/night with shrinkage) wins at month level; intervals calibrated on held-out residuals; NOAA-20 → S-NPP calibration.
- Space-time clustering of VIIRS detections into fire events; transparent rule-based types using ESA WorldCover and
  PEATMAP; Gaussian-mixture cross-check; manual NASA Worldview check.
- ERA5 100 m wind (Open-Meteo) to compute a daily upwind fire load for Kyiv (300 km, ±30°, distance decay, 0–2-day lags).
- Streamlit app reading only small processed parquet tables.

## Benefits
- Emergency responders and land managers see when and where burning is typically high (critical weeks), and how the
  current month compares with history — on a consistent scale across sensors.
- Scientists get a harmonized, documented, uncertainty-aware record with all validation numbers.
- Kyiv authorities get an evidence-based list of what burns upwind on poor-air days and where (stated as an association).

## Limitations
- Counts satellite detections, not officially reported fires; no detection ≠ no fire.
- Pre-2012 values are model estimates; small VIIRS-only fires cannot be placed in space before 2012.
- Burn types are model output; peat is under-detected near Kyiv; grass is hard to verify from imagery; manual check
  is small (18 scored events).
- The fire → PM link is not validated: CAMS misses the April 2020 smoke episode and is partly circular; ground
  stations exist only from Sept 2023.
- Conflict-affected areas are approximate oblast-level zones; trends after 2022 depend on war-time years; no causes
  are attributed; no forecasts.

## Tools
Python 3.13, pandas, pyarrow, NumPy, SciPy, statsmodels, scikit-learn, GeoPandas, Shapely, rasterio, Matplotlib,
Plotly, Streamlit, streamlit-folium / Folium; NASA FIRMS, NASA Worldview, ESA WorldCover, PEATMAP, Open-Meteo
(ERA5, CAMS), OpenAQ, geoBoundaries, OpenStreetMap; GitHub, Streamlit Community Cloud.

## Data sources
NASA FIRMS (MODIS C6.1 MCD14ML/MCD14DL; VIIRS 375 m VNP14IMG, VJ114IMG; NOAA-21) for Ukraine and Belarus · ESA
WorldCover 2021 v200 · PEATMAP (Xu et al. 2018, doi:10.5518/252, CC BY 4.0) · ERA5 via Open-Meteo · CAMS European
air quality via Open-Meteo · OpenAQ · geoBoundaries · OpenStreetMap (Chornobyl Exclusion Zone boundary).

## Use of AI tools
The code, analysis pipeline, app and documentation were developed with an AI coding assistant (Anthropic Claude,
via Claude Code) under the team's direction. The team defined the research questions, scientific rules and
validation requirements, made the method decisions (e.g. constant-area trends, conflict-zone handling), reviewed
the results, and manually labelled the 40 validation events in NASA Worldview. All numbers come from real data
processed by the code in this repository; no synthetic data were used.
