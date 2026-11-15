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
