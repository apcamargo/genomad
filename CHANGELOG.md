# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Dockerfile for version 1.0.0.
- `Sequence` class: add support for `str` in `__eq__`.
- `Sequence` class: add a `__hash__` method.

### Changed
- `Sequence` class: simplify `has_dtr` return statement.
- `Sequence` class: make `__repr__` more friendly for long sequences.
- `Sequence` class: rename the `id` property to `accession`.
- Amino acids are now written to `_provirus_aragorn.tsv`.

### Fixed
- tRNA coordinates are now 1-indexed.

## [1.0.0] - 2022-07-12
### Added
- First release.
