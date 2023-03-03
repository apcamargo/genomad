# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.5.0] - 2023-03-02
### Changed
- Given that geNomad applies a minimum score filter (since version 1.4.0), the help dialogue of the `--min-score` parameter was modified to remove the following sentence: *"By default, the sequence is classified as virus/plasmid if its virus/plasmid score is higher than its chromosome score, regardless of the value"*.
- The following parameters were added to the MMseqs2 search command: `--max-seqs 1000000 --min-ungapped-score 20 --max-rejected 225`. As a result, changing `--splits` won't affect the search results anymore.

## [1.4.0] - 2023-02-17
### Added
- Mention Docker and the NMDC EDGE implementation in the `README.md`.
- Add the `--min-plasmid-hallmarks-short-seqs` and `--min-virus-hallmarks-short-seqs` parameters. These options allow filtering out short sequences (less than 2,500 bp) that don't encode a minimum number of hallmark genes. By default, short sequences need to encode at least one hallmark to be classified as a virus or a plasmid.
- Add the `--conservative` and `--relaxed` presets that control post-classification filters. The `--conservative` option makes those filters even more aggressive, resulting in more restricted sets of plasmid and virus, containing only sequences whose classification is strongly supported. The `--relaxed` preset disables all post-classification filters.

### Changed
- Windows with more than 4,000 Ns are ignored when encoding sequences for the neural network classification. The first window is always processed, regardless of the amount of Ns.
- Changed the default value of `--min-score` from 0.0 to 0.7.
- Changed the default search sensitivity from 4.0 to 4.2.
- Update `README.md` to version 1.4.0. This includes mentions to the `--conservative` and `--relaxed` flags and a warning about how changes in `--splits` can affect geNomad's output.

## [1.3.3] - 2023-01-05
### Fixed
- Fix a bug in `score-calibration` that happened when `find-proviruses` was executed but no provirus was detected. The module now checks if proviruses were detected (using `utils.check_provirus_execution`) before counting the total number of sequences.

## [1.3.2] - 2022-12-28
### Fixed
- Require `numpy <1.6`. Fixes [#7](https://github.com/apcamargo/genomad/issues/7), which occurs because `numba` doesn't support `numpy >=1.24` yet.

## [1.3.1] - 2022-12-22
### Fixed
- Check if `find-proviruses` was executed when counting the number of sequences in the `score-calibration` module.

## [1.3.0] - 2022-12-12
### Added
- Add support for AMR annotation.
- Update database parsing to allow BUSCO-based USCGs.

### Changed
- Sequences with no terminal repeats will be flagged with `No terminal repeats`, as `Linear` can be misleading.
- Print the number of plasmids and viruses in the summary module.
- Set `click.rich_click.MAX_WIDTH` to `None`.
- Reduce the default `--sensitivity` to `4.0`.
- Update `README.md` to version 1.3.0.

### Fixed
- Set `prog_name` in `click.version_option`.

## [1.2.0] - 2022-11-15
### Changed
- Mention the Zenodo upload of geNomad's database in `README.md`.
- Add the following sentence for the help dialogue of the `--min-plasmid-marker-enrichment`, `--min-virus-marker-enrichment`, `--min-plasmid-hallmarks`, and `--min-virus-hallmarks` parameters: "This option will be ignored if the annotation module was not executed".
- Apply a uniform prior to the empirical sample composition in `score_batch_correction`. This will shrink the effect of calibration when the empirical composition distribution is very skewed.
- Reduce the `--min-score` in the `README.md` example to 0.7.

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
