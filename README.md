# Spatial prediction-powered inference: public code candidate

This code supports **Conditional validity of prediction-powered inference for spatially dependent finite populations**. This is the public code candidate; the separate data/reference bundle is not yet published. The research software is licensed under BSD-3-Clause (LICENSE); its scope is specified in LICENSE_SCOPE.md. Third-party data retain their own terms. Code repository: https://github.com/ShinShirota/spatial-ppi . A persistent data archive and DOI have not yet been assigned. See DATA_SOURCES.md and third_party_rights.md.

## Current availability

The 44-file code/documentation candidate is available in this repository. The data/reference ZIP and its SHA256SUMS.txt are available from [release v0.1.0-candidate](https://github.com/ShinShirota/spatial-ppi/releases/tag/v0.1.0-candidate). Use the release assets for the data and the tagged source archive for code. Full LUCAS replay and the combined manifest check require that separate bundle. The release tag pins the associated code; no DOI is claimed.

## Archive assembly

Extract `spatial-ppi-code-0.1.0-candidate.zip` and `spatial-ppi-data-0.1.0-candidate.zip` into the same parent directory. Both use the root `spatial-ppi-0.1.0-candidate/`; their files are disjoint. The data archive adds inputs/ and reference/, including a copy of the source notices in reference/. Run the commands below from that root. The code archive alone cannot reproduce the LUCAS analysis. Neither archive contains raw CSVs, TIFFs, manuscript files, virtual environments or development history.

Check each downloaded archive against SHA256SUMS.txt before extracting. After assembly, `python code/verify_manifest.py` verifies the combined package. The data archive is published as a release asset; code can be obtained from the release tag. Merge the extracted inputs/ and reference/ directories into the code root.

## Environment

Tested locally with Python 3.12.14 on macOS arm64 and Apple Accelerate; exact package versions are in requirements.txt and environment.json. Create a fresh virtual environment and install `python -m pip install -r requirements.txt`. This does not constitute third-party clean-environment validation.

Run from the package root. `inputs/` and `reference/` are a separate local data/reference bundle, intentionally excluded from normal Git tracking. The scripts resolve these relative to their own location; no author storage or account is needed. Download that bundle from the release and follow its source-data notices. Some directional outputs were retained as summaries only; the limitations below still apply.

## Verify and replay

Use a new output directory for every attempt. Prefix scientific commands with `python code/run_logged.py <new-log-directory>` to preserve start/end UTC, exit status, stdout and stderr. Keep failures and all attempted label sets. Do not relax tolerances after inspecting a mismatch.

```sh
python code/verify_manifest.py
python code/run_logged.py outputs/smoke-log python code/experiment.py --dataset lucas --replicates 10 --output outputs/smoke
python code/run_logged.py outputs/lucas-log python code/experiment.py --dataset lucas --output outputs/full_lucas
python code/summarize.py outputs/full_lucas
python code/compare.py reference/lucas outputs/full_lucas outputs/lucas_comparison.json
python code/run_logged.py outputs/mechanism-log python code/experiment.py --dataset mechanism --output outputs/mechanism
python code/summarize.py outputs/mechanism
python code/compare.py reference/mechanism outputs/mechanism outputs/mechanism_comparison.json
python code/run_logged.py outputs/gap-log python code/gap_run.py --output outputs/full_gap
python code/gap_summarize.py outputs/full_gap
python code/gap_validate.py outputs/full_gap
python code/compare_gap.py reference/gap outputs/full_gap outputs/gap_comparison.json
python code/run_logged.py outputs/directional-log python code/directional.py --fields 128 --labels 128 --seed 2026100901 --with-null --out-dir outputs/full_directional
python code/crosscondition.py outputs/full_lucas outputs/crosscondition
python code/build_tables.py --lucas outputs/full_lucas --mechanism outputs/mechanism --output outputs/tables_final
python code/audit_regression.py outputs/full_lucas outputs/audit_regression
python code/theory_bridge.py outputs/theory_bridge
python code/slope_boundary.py outputs/slope_boundary
```

For gap smoke testing use `gap_run.py --pilot --output <new-directory>`; for directional smoke use `directional.py --fields 2 --labels 2 --out-dir <new-directory>`. Smoke outputs never replace full outputs.

`check_tables.py outputs <manuscript-directory>` generates S4–S7 and compares all ten tables with main.tex and supplement.tex, including signed displayed numbers. The manuscript itself is not distributed in this code candidate. Frozen table excerpts are in reference/tables. `compare_additional.py outputs` checks directional, cross-condition and raw-frame reconstruction outputs (expects frame_rebuilt/ as below).

## Rebuild the analysis frame from public sources

Obtain the named LUCAS CSV from the official link in DATA_SOURCES.md. From the official WorldCover data-access page obtain **2021 v200 Map** tiles N57E021, N57E024 and N57E027. Do not use the 2020 product or a display image. Retain raw files separately; do not commit them. The local replay reused the existing source files, with no new download.

```sh
python code/rebuild_frame.py --csv raw/EE_2018_20200213.CSV.csv --tiles raw/ESA_WorldCover_10m_2021_v200_N57E021_Map.tif raw/ESA_WorldCover_10m_2021_v200_N57E024_Map.tif raw/ESA_WorldCover_10m_2021_v200_N57E027_Map.tif --output outputs/frame_rebuilt
python code/compare_additional.py outputs
```

See INPUT_DICTIONARY.md for response, prediction, coordinates, QR conventions and masking. Frozen label masks and audit references are supplied in the separate bundle; rebuilding the frame alone does not recreate those historical audit inputs.

## Interpretation and limitations

The fixed-frame target is a three-coefficient least-squares regression. Label sampling is conditional on that frame. Same-seed replay is distinct from independent SVD/derivative checks, and neither proves an asymptotic theorem. Preserve crossed or nested field/label dependence when interpreting MCSEs. True propensities/full outcomes are oracle diagnostics, not feasible tuning inputs. Signed slope weights are valid.

The directional script retains field summaries and aggregate controls, not individual label-set records. It reports all failed evaluations; this replay reported zero. The original limitation remains explicit. Numerical comparisons use the tolerances in specs/replay_spec_v1.json. Publication gates are listed in publication_gate.md.
