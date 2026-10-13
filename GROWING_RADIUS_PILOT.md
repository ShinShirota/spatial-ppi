The standalone archive and checksum are attached to [v0.2.0-candidate](https://github.com/ShinShirota/spatial-ppi/releases/tag/v0.2.0-candidate). Scientific scripts and saved outputs are in that asset; the automatically generated source ZIP alone does not include the pilot archive.

# Growing-radius Gaussian pilot: Supplement S6.6 / Table S8

Release v0.2.0-candidate adds this self-contained archive to the spatial PPI reproducibility package. The archive contains original scientific scripts, the frozen specification, 48 synthetic populations, all 3,072 realized label sets, all 9,216 method evaluations, and field-level/aggregate summaries. It requires no LUCAS data, external rasters, private account, or manuscript TeX.

## Verify the archived table (Python standard library only)

Extract spatial-ppi-growing-pilot-0.2.0-candidate.zip and enter its directory:

```sh
python verify.py
```

This checks the distribution SHA256 manifest and reaggregates all Table S8 coverage/MCSE cells. It does not generate new samples. MCSE is the sample standard deviation of 16 population-specific coverage rates divided by sqrt(16), in percentage points. Failed and undefined evaluations stay in the denominator; there were zero of each among 9,216 attempts.

## Validate or reproduce

Install requirements.txt in a separate virtual environment. From the extracted directory use a NEW output directory for every action:

```sh
python launch.py validate outputs/algebra-check
python launch.py reaggregate outputs/saved-table-check
python launch.py replay outputs/full-replay
```

Validation uses the archived finite-algebra/derivative diagnostic seeds and tolerances. Replay executes validation first, then all 48 populations x 64 label acquisitions x 3 methods, using master seed 2026101302. It can require about 0.6 GiB free disk and 1 GiB RAM; it preserves logs, failures, and every attempted sample. The launcher copies the unchanged scripts into a new working tree with src/regression_core.py; historical results are never overwritten. Same-seed replay is not independent implementation validation. Cross-platform floating-point equality is not guaranteed; the archived algebra and finite-difference tolerances are unchanged. Do not silently relax them. Publication packaging checked the finite-algebra validation and saved-output reconciliation, not a fresh full scientific replay.

## Table mapping and limits

| Manuscript | Source | Generation / check |
|---|---|---|
| S6.6 / Table S8 | pilot/results/replicates.jsonl, field_summary.json, summary.json | pilot/run.py; pilot/summarize.py; verify.py |
| S6.6 gap/remainder/covariance diagnostics | same saved per-evaluation records and summaries | pilot/engine.py and pilot/summarize.py |

At n=128,512,2048, hard slope coverage is 90.72,93.16,93.55%; soft+BC coverage is 88.77,92.87,93.75%; known-pi HT coverage is 96.00,95.51,94.82%. The targeted methods do not attain 95% at the largest level. Levels have independent populations and changing acquisition modulation/structural variance. This exploratory pilot is not a proof of conditional asymptotic coverage, a test of the directional-versus-energy dimensional gain, or a non-Gaussian validation.

## Provenance, rights and chronology

Scientific scripts, spec_v1.json, validation.json, numerical JSON/JSONL and NPZ files are retained unchanged, except for distribution metadata described next. Private machine directory prefixes in distributed provenance input_hashes/environment and execution logs are replaced by <PROJECT_ROOT>. Original historical output_hashes.json is preserved as an audit record; use the top-level manifest.json for distribution integrity because normalized metadata has different bytes. All original files remain retained by the author.

The original spec and Japanese report describe the earlier, unpublished pilot stage. Those historical statements are not current publication restrictions: the author subsequently adopted S6.6/Table S8 and authorized this public release on 2026-10-13. No seeds, outcomes, failure records or scientific tolerances were changed for publication. The current distribution statement is this README.

Code, documentation and author-generated synthetic pilot outputs are provided under the included BSD-3-Clause license. The old LUCAS/WorldCover data retain their original third-party terms and are not duplicated here. No raw third-party data, manuscript, reviewer messages, credentials or unrelated local files are included. Earlier release tags and assets remain unchanged. No DOI is assigned.
