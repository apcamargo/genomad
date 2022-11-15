# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- Mention the Zenodo upload of geNomad's database in `README.md`.
- Add the following sentence for the help dialogue of the `--min-plasmid-marker-enrichment`, `--min-virus-marker-enrichment`, `--min-plasmid-hallmarks`, and `--min-virus-hallmarks` parameters: "This option will be ignored if the annotation module was not executed".
- Apply a uniform prior to the empirical sample composition in `score_batch_correction`. This will shrink the effect of calibration when the empirical composition distribution is very skewed.

### Fixed
- Fix a bug in the score calibration module where the sample size was set to a constant value and the "Your sample has less than 1,000 sequences…" warning would always appear.

## [1.1.0] - 2022-08-22
### Added
- Dockerfile for version 1.0.0.
- `Sequence` class: add support for `str` in `__eq__`.
- `Sequence` class: add a `__hash__` method.
- Compute marker enrichment in the `marker-classification` module.
- Add columns for plasmid and virus marker enrichment to the `_plasmid_summary.tsv` and `_virus_summary.tsv` files.
- Set `--min-plasmid-marker-enrichment` and `--min-virus-marker-enrichment` to `0` as default. This will alter the results when using default parameters.
- Add support for plasmid and virus hallmarks. Requires geNomad database v1.1.
- Add CONJscan annotations to `_plasmid_summary.tsv`. Requires geNomad database v1.1.

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
- Fix a problem in `DatabaseDownloader.get_version` where only the major version was compared.

## [1.0.0] - 2022-07-12
### Added
- First release.
