# Inputs and fixed analysis frame

All 2,665 Estonia LUCAS 2018 records are kept in original order. POINT_ID is a public site identifier, not anonymized. Coordinates identify public sample locations; the inspected bundle is eligible for redistribution under the recorded source terms, pending author authorization for publication. Coordinates derive from public theoretical-point TH_LONG/TH_LAT, not surveyor GPS or parking positions.

- `inputs/inputs/lucas_worldcover_joined.csv`: POINT_ID; easting_3035, northing_3035 (metres, EPSG:3035).
- `inputs/source_snapshot/spatial_ppi_stage11/results/primary/population_inputs.npz`: Y = woodland indicator (LC1 begins with C); f = WorldCover 2021 v200 native Map code 10; B = 15-column polynomial QR basis; Bmap = 16-column map/polynomial QR basis; POINT_ID = matching site order. Rebuild code is authoritative for column ordering and positive-diagonal QR convention.
- The nine other NPZ inputs preserve original Bernoulli label masks, sizes and true inclusion probabilities, plus old intercept estimates/variances for the numerical audit. Probabilities and full-population outcomes are oracle inputs only where explicitly designated by the scientific scripts.
- Scientific X is always [1, centred easting/100000, centred northing/100000]. It is unchanged across methods. B and Bmap are calibration dictionaries, distinct from X.
- `inputs/protocol/protocol.json` and the stage11 protocol preserve original budgets, tuning and sampling specifications. Frozen label masks define the exact replay; no samples are redrawn.
- No height rasters, canopy data, photographs or neighbourhood features are needed for these analyses.

`rebuild_frame.py` checks unique IDs, exactly one covering tile, missing Map pixels and QR rank. A mismatch is an error, not an exclusion rule. Original CSV footer text is distinguished from records. All scientific arrays and coordinates rebuilt from the held raw files passed the predeclared comparison tolerances.
