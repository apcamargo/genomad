# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Dockerfile for version 1.0.0.
- `Sequence` class: add support for `str` in `__eq__`.
- `Sequence` class: add a `__hash__` method.
- Compute marker enrichment in the `marker-classification` module.
- Add columns for plasmid and virus marker enrichment to the `_plasmid_summary.tsv` and `_virus_summary.tsv` files.
- Set `--min-plasmid-marker-enrichment` and `--min-virus-marker-enrichment` to `0` as default. This will alter the results when using default parameters.
- Add support for plasmid and virus hallmarks. Requires geNomad datavase v1.1.

### Changed
- `Sequence` class: simplify `has_dtr` return statement.
- `Sequence` class: make `__repr__` more friendly for long sequences.
- `Sequence` class: rename the `id` property to `accession`.
- Amino acids are now written to `_provirus_aragorn.tsv`.
- Update the XGBoost model file to the `.ubj` format.
- Require `xgboost >=1.6`.
- The taxonomic lineage in `_taxonomy.tsv` and `_virus_summary.tsv` will use `Viruses` as the highest rank, instead of `root`.
- Change order of the columns in `_plasmid_summary.tsv` and `_virus_summary.tsv`.
- Explicitly set `fraction` to `0.5` in `taxopy.find_majority_vote`.

### Fixed
- tRNA coordinates are now 1-indexed.
- Write `summary_execution_info`.
- Fix a problem in `DatabaseDownloader.get_version` where it only only compared the major version.

## [1.0.0] - 2022-07-12
### Added
- First release.
