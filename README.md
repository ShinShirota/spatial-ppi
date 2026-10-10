# Spatial prediction-powered inference: reproducibility package

Code for **Conditional validity of prediction-powered inference for spatially dependent finite populations**. Version 0.1.1-candidate corrects distribution documentation only; scientific code and the data bundle are unchanged from v0.1.0-candidate. The software uses BSD-3-Clause; third-party data retain their own terms. No DOI is assigned.

## Current availability and download

Download **Source code (zip)** from this version's GitHub Release, or check out its tag. The v0.1.1-candidate source ZIP normally expands to `spatial-ppi-0.1.1-candidate/`.

The unchanged data ZIP and its separate checksum asset are available from [the v0.1.0-candidate data release](https://github.com/ShinShirota/spatial-ppi/releases/tag/v0.1.0-candidate):

- `spatial-ppi-data-0.1.0-candidate.zip`
- `SHA256SUMS.txt` (a separate release asset, not inside either ZIP)

SHA256SUMS.txt covers **the data ZIP only**, not GitHub's automatically generated source archives. The data ZIP SHA-256 is `78a6f83dcc4f876676ad49a759c713f7881a33b76c84a8d64959fc8da2e91f5d` (139,665,116 bytes). These hashes check integrity; they are not an external signature.

## Assembly

Place both ZIPs and SHA256SUMS.txt in a new directory. On macOS:

```sh
shasum -a 256 -c SHA256SUMS.txt
unzip spatial-ppi-0.1.1-candidate.zip
unzip spatial-ppi-data-0.1.0-candidate.zip
cp -R spatial-ppi-0.1.0-candidate/inputs spatial-ppi-0.1.1-candidate/
cp -R spatial-ppi-0.1.0-candidate/reference spatial-ppi-0.1.1-candidate/
cd spatial-ppi-0.1.1-candidate
python code/verify_manifest.py
```

GitHub may name the downloaded source ZIP differently; use its actual filename. The data archive retains its original version/root intentionally. The combined manifest verifies 113 recorded files; it is not an external signature. Source rasters, photographs, manuscript files and virtual environments are excluded. Preserve the source-rights notices in `reference/`.

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

## Table checks without manuscript TeX

`check_tables.py` compares generated tables against the author's manuscript TeX, which is not distributed. A user without that TeX can instead compare the seven generated table files against the frozen reference files after the full LUCAS/mechanism runs and `build_tables.py` above:

```sh
python - <<'CHECK'
from pathlib import Path
names = ['mechanism', 'lucas_east', 'joint', 'lucas_north', 'paired', 'smalln', 'truth']
failed = []
for name in names:
    filename = name + '_table.tex'
    ok = (Path('outputs/tables_final') / filename).read_bytes() == (Path('reference/tables') / filename).read_bytes()
    print(filename, 'PASS' if ok else 'FAIL')
    if not ok:
        failed.append(filename)
raise SystemExit(bool(failed))
CHECK
```

This checks Main Tables 1–3, Supplement Tables S1–S3, and the auxiliary truth table against the distributed reference; it does not independently validate their transcription into the manuscript. S4–S6 frozen display rows are also in `reference/tables/`. S7 is represented by the machine-readable cross-condition summaries in `reference/crosscondition/`; there is no standalone frozen S7 TeX file. The gap and additional-output comparisons above check their numerical sources. They should not be described as an automatic comparison of all ten tables to an undistributed manuscript.
