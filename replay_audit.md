# Local replay audit — 2026-10-10

All ten manuscript tables (119 rows) match regenerated values, including original signs and displayed digits. Main1–3/S1–3 were regenerated from full estimator outputs; S4 from the full fixed-dictionary experiment; S5/S6 from the full directional replay; S7 from reaggregating the new full LUCAS output. A first formatting comparison identified only an old plain-X method label; the table generator now uses the manuscript's bold X. No numerical values or manuscript files were changed.

- LUCAS: 199 array/summary-column comparisons pass. All 18,000 saved label sets; 216,000 method evaluations, both prediction variants. Fit exceptions: zero; undefined joint regions remain separate from fit exceptions.
- Constructed mechanism: 161 comparisons pass; 32 fields, 500 shared label draws in each original cell.
- Fixed-dictionary gap: all 36,864 attempts and 917,733 float values compared; maximum absolute difference zero. Original failure/non-reporting denominators retained.
- Directional: 128 fields × 128 labels, 589,824 primary and 196,608 control evaluations; zero reported failures. CSV comparison normalizes method names and column order, preserving values and predeclared rounding tolerance. Individual label-level records were not saved by the original script.
- Additional checks: 111 column/array comparisons pass across directional results, both LUCAS cross-condition variants, and rebuilt source frame.
- Independent SVD: 297 cases; maximum estimate error 4.774e-15, covariance error 6.245e-17. 33 finite differences, 33 deletions, 60 exact-remainder examples and 1,242 summary rows checked; audit reports no failures.
- Deterministic theory/code and slope-boundary identities pass their original tolerances. These finite checks are not asymptotic proofs.

LUCAS, gap and directional smoke tests were separate from the full runs. Initial directional smoke failed before simulation because threadpoolctl was absent. Initial dependency download failed under network restrictions; an approval attempt timed out; the permitted retry installed threadpoolctl 3.6.0 in the project environment and the new smoke attempt passed. All failed logs are preserved internally. An initial environment probe also failed on missing package metadata and was corrected without altering results.

Actual UTC start/end times, commands, process and exit statuses are retained in the internal execution index; source/code/input hashes and original outputs are retained in the local audit directory. Public environment versions are in environment.json and requirements.txt. This is local same-seed replay plus specified independent numerical checks, not third-party clean-environment reproduction.
