> Current distribution note (10 October 2026): the inspected data bundle is published with author authorization in GitHub release v0.1.0-candidate. Any pending-publication wording below records the earlier rights review, not the current release status. Source conditions are unchanged.

# Third-party rights and release decision — 2026-10-10

Official sources checked on 10 October 2026:

- Eurostat: https://ec.europa.eu/eurostat/help/copyright-notice . Statistical reuse requires source acknowledgement; individual notices and third-party exceptions apply. Modifications require disclosure and a Eurostat non-responsibility statement.
- WorldCover: https://esa-worldcover.org/en/data-access ; dataset DOI https://doi.org/10.5281/zenodo.7254221 . Retain the product's CC BY 4.0 attribution and distinguish extracted/modified data from the source product.

Author-code licensing does not replace those terms. The research software now uses BSD-3-Clause; LICENSE_SCOPE.md excludes third-party data from that software grant. Copyright ownership and permission granted to downstream users are separate questions.

**Decision after file-level review:** the existing minimal LUCAS/WorldCover analysis and reference bundle is eligible for redistribution with the attribution, modification notices and source terms recorded here and in DATA_SOURCES.md. This is a reasoned assessment of the published terms and inspected files, not individual clearance from Eurostat or ESA. The previous unresolved data-rights gate is closed for this exact candidate content; author authorization for actual publication remains required.

The official Estonia CSV downloaded on 10 October 2026 is byte-identical to the held input (SHA256 21e25a725b63ad36c5b735c24cd8527af2dbff4a8c743dd73dadd0470d0e5631). All 2,665 POINT_ID values and the woodland indicator match; EPSG:3035 coordinates are transformations of public TH_LONG/TH_LAT, agreeing within 9.78e-9 metres. The official WEB CSV Record Descriptor identifies these as theoretical point coordinates, separately from survey GPS and parking coordinates. No special restriction on these fields was found in the reviewed official distribution page, descriptor or reuse notice. The separately requested grassland dataset is not used.

Scope: `inputs/inputs/lucas_worldcover_joined.csv`, the numerical arrays in `inputs/source_snapshot/`, and LUCAS-derived `reference/lucas/`, `reference/crosscondition/` and table excerpts. Sampling masks, estimator outputs, specifications and simulated results were inspected by file/schema; the candidate does not add owner names, addresses, photographs, parking/GPS observation fields, restricted records or private annotations. Site IDs and locations remain identifiable public geographic data, not anonymized data. This finding does not license future personal-data linkage or other datasets.

WorldCover-derived tree indicators retain CC BY 4.0 source attribution, the licence link and modification notice. CC BY 4.0 permits sharing adaptations subject to its conditions; it does not transfer trademark, privacy or other excluded rights. Do not apply the research-code licence as a replacement for either source's terms. Retain these notices with the data bundle. Raw CSV and source rasters remain excluded by packaging choice, not because redistribution of all source data is categorically forbidden.

Additional evidence:
- https://ec.europa.eu/eurostat/documents/205002/8072634/LUCAS2018-RecordDescriptor-190611.pdf
- https://creativecommons.org/licenses/by/4.0/legalcode.en


Credit: Eurostat, LUCAS 2018, Estonia. Derived indicators, coordinate transformations and analyses are the author's responsibility; Eurostat is not responsible for the modifications or interpretation. ESA WorldCover project 2021; Copernicus Sentinel data processed by the ESA WorldCover consortium; Zanaga et al. (2022), ESA WorldCover 10 m 2021 v200.
