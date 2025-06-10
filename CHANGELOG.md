# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.11.1] - 2025-06-10
### Changed
- Simplify the `softmax` function implementation and add an `axis` parameter.
- Use `os.sched_getaffinity` instead of `multiprocessing.cpu_count()` to determine the number of available CPUs.

## [1.11.0] - 2025-02-18
### Added
- Add the `--lenient-taxonomy` and `--full-ictv-lineage` options to the `annotate`, `find-proviruses`, and `end-to-end` modules. The `--lenient-taxonomy` option enables assignment of genomes to taxa below the family rank (subfamily, genus, subgenus, and species). The `--full-ictv-lineage` option enables the output of the full ICTV lineage of each genome, including ranks that are hidden by default (subrealm, subkingdom, subphylum, subclass, suborder, subfamily, and subgenus).

### Changed
- Remove the `--conservative-taxonomy` option of the `annotate` and `end-to-end` modules.

## [1.10.0] - 2025-02-15
### Changed
- Add support for the geNomad database v1.8, which introduced an additional column in the metadata table.
- Update `README.md` to the database v1.8.
- Set the minimum Python version to `3.9`.

## [1.9.0] - 2025-01-25
### Changed
- Add installation instructions using Pixi.
- Use raw strings for the regular expression in `utils.natsort`.
- Update the version requirements for `keras` to `>=3` and `tensorflow` to `>=2.16`.

## [1.8.1] - 2024-10-05
### Fixed
- Write the `min_number_genes` value to the parameters JSON file of the `summary` module.
- Set maximum `tensorflow` version to below `2.16`.

### Changed
- Set the `break_on_hyphens` parameter of the `textwrap.fill` function to `False` to prevent line breaks at `-` characters. This ensures that sequences with gaps in FASTA files generated using `Sequence.__str__()` maintain consistent line width.
- Compare `Enum` by identity in the `open_file` function.

## [1.8.0] - 2024-04-10
### Added
- Added the `--min-number-genes` parameter to the `summary` module. This parameter allows users to set the minimum number of genes a sequence must encode to be considered for classification as a plasmid or virus. The default value is `1`. When `--conservative` is used, this parameter is set to `1`. When `--relaxed` is used, this parameter is set to `0`. This filter has no effect if the `annotate` module is not executed.

### Changed
- Added a hyperlink to the official documentation in the help dialogue.
- The virus taxonomic lineage is presented using a fixed number of fields separated by semicolons (`;`). As a result, for genomes that could not be assigned to the family level (the most specific taxonomic rank), there will be trailing semicolons at the end of the lineage string.
- Do not apply the gene-based post-classification filters when the `annotate` module is not executed.
- Set the default value of `--min-plasmid-marker-enrichment` to `0.1`.

## [1.7.6] - 2024-03-19
### Fixed
- Set maximum `keras` version to below `3.0`. This prevents errors due to incompatibility with `keras >=3.0`, such as the `shape` parameter not accepting an integer as input.

## [1.7.5] - 2024-03-03
### Fixed
- Set the `CUDA_VISIBLE_DEVICES` environment variable to `-1` in `nn_classification`. This fixes a bug where the `nn_classification` module would fail to run when a GPU was available and the input had a single sequence.

## [1.7.4] - 2023-12-08
### Fixed
- Fixed the parsing of MMseqs2 integrase output to extract only the gene accession, rather than the entire header. This addresses a bug introduced in version 1.5.2, where the integrase gene accession was not accurately parsed because the entire header was extracted. As a result, the `find-proviruses` module can now properly add integrases to gene tables and extend boundaries using integrase coordinates.

### Changed
- Replace ambiguous variable name in `read_fasta`.
- Define name `current_contig` at the beginning of `_append_aragorn_tsv`.

## [1.7.3] - 2023-11-30
### Fixed
- Set minimum `pyrodigal-gv` version to `0.3.1`. This fixes a bug introduced in `0.3.0` that led to the identification of RBS motifs not reported by Prodigal.

### Changed
- Remove the `CCGGGG` RBS motif from the list of motifs.

## [1.7.2] - 2023-11-28
### Fixed
- Add the `CCGGGG` RBS motif to the list of motifs.

### Changed
- Do not include stop codon (`*`) at the end of protein sequences.
- Set minimum `pyrodigal-gv` version to `0.2.0`.

## [1.7.1] - 2023-10-26
### Changed
- Replace `prodigal-gv` with `pyrodigal-gv`

## [1.7.0] - 2023-09-13
### Changed
- The `mmseqs search` command has been replaced by a two-step alignment workflow. In the first alignment step, `--alignment-mode 1` and `--max-rejected` are utilized, while the second step uses `--alignment-mode 2` and `-c 0.2`. This change reduces the number of alignments that are rejected due to not meeting the minimum coverage cutoff and mitigates the issue where the annotation results change when the input sequence order is altered.
- The `--min-ungapped-score` parameter of `mmseqs prefilter` was increased from `20` to `25`.
- The `--max-rejected` parameter of the first `mmseqs align` step was increased from `225` to `280`.

## [1.6.1] - 2023-07-31
### Fixed
- Replace `np.warnings` with `warnings` to add compatibility with `numpy >= 1.24`.

## [1.6.0] - 2023-07-31
### Changed
- Update `numba` (`>=0.57`) and `numpy` (`>=1.21`) version requirements.
- Use `casefold` for sequence comparison within the `Sequence` class.
- Remove type annotations of methods of the `Sequence` class that return an instance of `Sequence`.
- Use `console.status` to log the deletion of the `.tar.gz` file during the execution of `download-database`.
- Make the conservative assignment at the family level optional via the `--conservative-taxonomy` parameter. This increases the amount of viral genomes assigned to a family when executing geNomad with default parameters.

### Fixed
- Fix parameter names in the error message of `--conservative` and `--relaxed` (e.g. `--min_score` → `--min-score`).

## [1.5.2] - 2023-05-11
### Added
- Display a progress bar showing the progress of the classification process in `nn-classification`.

### Changed
- Update `README.md` to the database version 1.3.0.

### Fixed
- Make `mmseqs convertalis` output the whole sequence header instead of gene accesions. This prevents parsing conflits with geNomad's other components in cases where MMseqs2 uses its built-in special parsers for specific header formats (e.g. RefSeq).

## [1.5.1] - 2023-03-30
### Added
- Add the `--threads` parameter to the `nn-classification` module, which allows controlling the number of threads used for classifying sequences using the neural network model.

### Changed
- Mention post-classification filters the in the `summary` module description.

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
