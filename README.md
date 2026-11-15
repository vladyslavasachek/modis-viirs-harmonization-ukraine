# Harmonized burning-activity calendar for Ukraine (MODIS + VIIRS)

NASA Space Apps Challenge 2026 — *Harmonization of MODIS and VIIRS Hot Spots*.
Full README (problem, data, method, validation, results, limitations, how to run, attribution) is completed in
Stage 8; see `PROGRESS.md` for the full log.

## Burn-type validation (manual Worldview check)

40 fire events (8 per predicted type, stratified random among VIIRS S-NPP events with ≥ 3 detections, 2012–2025,
outside conflict-affected zones) were labelled by one team member in NASA Worldview (true colour during the fire
and burn scar after it); the file is `validation_events_worldview.csv`.

- **22 of 40** events could be labelled confidently; **18 were "unclear"** from the imagery.
- On the **18** confident labels whose predicted type is a definite class, **14 agree**:
  **accuracy 78% (95% Wilson interval 55–91%)**. The `other / uncertain` class makes no type claim; its 4 confident
  labels were 3 field and 1 forest.
- By predicted type: crop residue 7/7 correct (1 unclear of 8), forest 4/4 (4 unclear), grass 1/1 (**7 unclear**),
  peat **2/6** (2 unclear) — all four disagreements are in the peat class (2 looked like fields, 2 like forest,
  possibly forested peat or field/forest edges).
- Adding PEATMAP did not change any of the 40 sampled labels.
- **Framing:** small, stratified (not prevalence-weighted), single-labeller sample → an indicative check, not an
  accuracy estimate for the whole map. Crop residue and forest labels look reliable; grass is essentially
  unvalidated from imagery; peat labels are the weakest class.
