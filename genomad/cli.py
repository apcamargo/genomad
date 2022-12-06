import multiprocessing
from pathlib import Path

import genomad
import rich_click as click


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.MAX_WIDTH = None
click.rich_click.STYLE_OPTIONS_TABLE_BOX = "SIMPLE"
click.rich_click.STYLE_COMMANDS_TABLE_SHOW_LINES = True
click.rich_click.STYLE_COMMANDS_TABLE_PAD_EDGE = True
click.rich_click.STYLE_COMMANDS_TABLE_BOX = "SIMPLE"

click.rich_click.COMMAND_GROUPS = {
    "genomad": [
        {
            "name": "Database download",
            "commands": [
                "download-database",
            ],
        },
        {
            "name": "End-to-end execution",
            "commands": [
                "end-to-end",
            ],
        },
        {
            "name": "Modules",
            "commands": [
                "annotate",
                "find-proviruses",
                "marker-classification",
                "nn-classification",
                "aggregated-classification",
                "score-calibration",
                "summary",
            ],
        },
    ]
}

click.rich_click.OPTION_GROUPS = {
    "genomad download-database": [
        {
            "name": "Basic options",
            "options": ["--keep", "--verbose"],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
    "genomad annotate": [
        {
            "name": "Basic options",
            "options": ["--cleanup", "--restart", "--threads", "--verbose"],
        },
        {
            "name": "Advanced options",
            "options": ["--sensitivity", "--evalue", "--splits", "--use-minimal-db"],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
    "genomad find-proviruses": [
        {
            "name": "Basic options",
            "options": [
                "--cleanup",
                "--restart",
                "--skip-integrase-identification",
                "--skip-trna-identification",
                "--threads",
                "--verbose",
            ],
        },
        {
            "name": "Advanced options",
            "options": [
                "--crf-threshold",
                "--marker-threshold",
                "--marker-threshold-integrase",
                "--marker-threshold-edge",
                "--max-integrase-distance",
                "--max-trna-distance",
                "--sensitivity",
                "--evalue",
            ],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
    "genomad marker-classification": [
        {
            "name": "Basic options",
            "options": ["--restart", "--threads", "--verbose"],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
    "genomad nn-classification": [
        {
            "name": "Basic options",
            "options": ["--cleanup", "--restart", "--verbose"],
        },
        {
            "name": "Advanced options",
            "options": ["--single-window", "--batch-size"],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
    "genomad aggregated-classification": [
        {
            "name": "Basic options",
            "options": ["--restart", "--verbose"],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
    "genomad score-calibration": [
        {
            "name": "Basic options",
            "options": ["--composition", "--verbose"],
        },
        {"name": "Advanced options", "options": ["--force-auto"]},
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
    "genomad summary": [
        {
            "name": "Basic options",
            "options": [
                "--min-score",
                "--max-fdr",
                "--min-plasmid-marker-enrichment",
                "--min-virus-marker-enrichment",
                "--min-plasmid-hallmarks",
                "--min-virus-hallmarks",
                "--max-uscg",
                "--verbose",
            ],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
    "genomad end-to-end": [
        {
            "name": "Basic options",
            "options": [
                "--cleanup",
                "--restart",
                "--disable-find-proviruses",
                "--disable-nn-classification",
                "--enable-score-calibration",
                "--threads",
                "--verbose",
            ],
        },
        {
            "name": "annotation options",
            "options": ["--sensitivity", "--splits"],
        },
        {
            "name": "find-proviruses options",
            "options": [
                "--skip-integrase-identification",
                "--skip-trna-identification",
            ],
        },
        {
            "name": "score-calibration options",
            "options": ["--composition", "--force-auto"],
        },
        {
            "name": "summary options",
            "options": [
                "--min-score",
                "--max-fdr",
                "--min-plasmid-marker-enrichment",
                "--min-virus-marker-enrichment",
                "--min-plasmid-hallmarks",
                "--min-virus-hallmarks",
                "--max-uscg",
            ],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(prog_name="geNomad")
def cli():
    """
    [cyan]geNomad[/cyan]: Identification of mobile genetic elements

    Read the documentation at: https://portal.nersc.gov/genomad/
    """


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("destination", type=click.Path(path_type=Path, exists=True))
@click.version_option(prog_name="geNomad")
@click.option(
    "--keep",
    is_flag=True,
    default=False,
    show_default=True,
    help="Do not delete the compressed database file.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
def download_database(destination, keep, verbose):
    """
    Download the latest version of geNomad's database and save it in the
    [u]DESTINATION[/u] directory.
    """
    genomad.download.main(destination, keep, verbose)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.argument("database", type=click.Path(path_type=Path, exists=True))
@click.version_option(prog_name="geNomad")
@click.option(
    "--cleanup",
    is_flag=True,
    default=False,
    show_default=True,
    help="Delete intermediate files after execution.",
)
@click.option(
    "--restart",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite existing intermediate files.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="Number of threads to use.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
@click.option(
    "--sensitivity",
    "-s",
    type=click.FloatRange(min=0.0),
    default=4.0,
    show_default=True,
    help="""MMseqs2 marker search sensitivity. Higher values will annotate more
            proteins, but the search will be slower and consume more memory.""",
)
@click.option(
    "--evalue",
    "-e",
    type=float,
    default=1e-3,
    show_default=True,
    help="Maximum accepted E-value in the MMseqs2 search.",
)
@click.option(
    "--splits",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="""Split the data for the MMseqs2 search. Higher values will reduce
            memory usage, but will make the search slower. If the MMseqs2 search
            is failing, try to increase the number of splits.""",
)
@click.option(
    "--use-minimal-db",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Use a smaller marker database to annotate proteins. This will make
            execution faster but sensitivity will be reduced.""",
)
def annotate(
    input,
    output,
    database,
    use_minimal_db,
    restart,
    threads,
    verbose,
    sensitivity,
    evalue,
    splits,
    cleanup,
):
    """
    Predict the genes in the [u]INPUT[/u] file (FASTA format), annotate them
    using geNomad's markers (located in the [u]DATABASE[/u] directory), and
    write the results to the [u]OUTPUT[/u] directory.
    """
    genomad.annotate.main(
        input,
        output,
        database,
        use_minimal_db,
        restart,
        threads,
        verbose,
        sensitivity,
        evalue,
        splits,
        cleanup,
    )


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.argument("database", type=click.Path(path_type=Path, exists=True))
@click.version_option(prog_name="geNomad")
@click.option(
    "--restart",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite existing intermediate files.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
@click.option(
    "--cleanup",
    is_flag=True,
    default=False,
    show_default=True,
    help="Delete intermediate files after execution.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="Number of threads to use.",
)
@click.option(
    "--skip-integrase-identification",
    is_flag=True,
    default=False,
    show_default=True,
    help="Disable provirus boundary extension using nearby integrases.",
)
@click.option(
    "--skip-trna-identification",
    is_flag=True,
    default=False,
    show_default=True,
    help="Disable provirus boundary extension using nearby tRNAs.",
)
@click.option(
    "--crf-threshold",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0.4,
    show_default=True,
    help="""Minimum gene-level score to flag a provirus gene using the conditional
            random field model. Lower values will result in longer proviruses but
            will increase the probability of host genes being flagged as part of
            proviruses.""",
)
@click.option(
    "--marker-threshold",
    type=float,
    default=12.0,
    show_default=True,
    help="""Minimum total virus marker score allowed for proviruses that do not
            encode integrases or are not located at scaffold edges. Lower
            values will increase the sensitivity but reduce the precision of
            the provirus identification procedure.""",
)
@click.option(
    "--marker-threshold-integrase",
    type=float,
    default=8,
    show_default=True,
    help="""Minimum total virus marker score allowed for proviruses that
            encode integrases.""",
)
@click.option(
    "--marker-threshold-edge",
    type=float,
    default=8,
    show_default=True,
    help="""Minimum total virus marker score allowed for proviruses that are
            located at scaffold edges.""",
)
@click.option(
    "--sensitivity",
    "-s",
    type=click.FloatRange(min=0.0),
    default=8.2,
    show_default=True,
    help="""MMseqs2 integrase search sensitivity. Higher values will identify more
            integrases, but the search will be slower and consume more memory.""",
)
@click.option(
    "--evalue",
    "-e",
    type=float,
    default=1e-3,
    show_default=True,
    help="Maximum accepted E-value in the MMseqs2 integrase search.",
)
@click.option(
    "--max-integrase-distance",
    type=click.IntRange(min=0),
    default=10_000,
    show_default=True,
    help="""Maximum allowed distance between provirus boundaries and the integrases
            used for boundary extension.""",
)
@click.option(
    "--max-trna-distance",
    type=click.IntRange(min=0),
    default=5_000,
    show_default=True,
    help="""Maximum allowed distance between provirus boundaries and the tRNAs used
            for boundary extension.""",
)
def find_proviruses(
    input,
    output,
    database,
    cleanup,
    restart,
    skip_integrase_identification,
    skip_trna_identification,
    threads,
    verbose,
    crf_threshold,
    marker_threshold,
    marker_threshold_integrase,
    marker_threshold_edge,
    max_integrase_distance,
    max_trna_distance,
    sensitivity,
    evalue,
):
    """
    Find integrated viruses within the sequences in [u]INPUT[/u] file
    using the geNomad markers (located in the [u]DATABASE[/u] directory) and
    write the results to the [u]OUTPUT[/u] directory. This command depends on
    the data generated by the [cyan]annotate[/cyan] module.
    """
    genomad.find_proviruses.main(
        input,
        output,
        database,
        cleanup,
        restart,
        skip_integrase_identification,
        skip_trna_identification,
        threads,
        verbose,
        crf_threshold,
        marker_threshold,
        marker_threshold_integrase,
        marker_threshold_edge,
        max_integrase_distance,
        max_trna_distance,
        sensitivity,
        evalue,
    )


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.argument("database", type=click.Path(path_type=Path, exists=True))
@click.version_option(prog_name="geNomad")
@click.option(
    "--restart",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite existing intermediate files.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="Number of threads to use.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
def marker_classification(input, output, database, restart, threads, verbose):
    """
    Classify the sequences in the [u]INPUT[/u] file (FASTA format) based on the
    presence of geNomad markers (located in the [u]DATABASE[/u] directory) and
    write the results to the [u]OUTPUT[/u] directory. This command depends on
    the data generated by the [cyan]annotate[/cyan] module.
    """
    genomad.marker_classification.main(
        input, output, database, restart, threads, verbose
    )


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.version_option(prog_name="geNomad")
@click.option(
    "--restart",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite existing intermediate files.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
@click.option(
    "--cleanup",
    is_flag=True,
    default=False,
    show_default=True,
    help="Delete intermediate files after execution.",
)
@click.option(
    "--single-window",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Use only the first window (6,000 bases) of each sequence to perform
            the classification. This will make execution faster and reduce memory
            usage, but the classification accuracy will decrease.""",
)
@click.option(
    "--batch-size",
    type=int,
    default=128,
    show_default=True,
    help="""Number of data points per batch of prediction. Use a smaller value
            to reduce memory comsumption at the cost of speed.""",
)
def nn_classification(
    input, output, single_window, batch_size, restart, verbose, cleanup
):
    """
    Classify the sequences in the [u]INPUT[/u] file (FASTA format) using the
    geNomad neural network and write the results to the [u]OUTPUT[/u] directory.
    """
    genomad.nn_classification.main(
        input, output, single_window, batch_size, restart, verbose, cleanup
    )


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.version_option(prog_name="geNomad")
@click.option(
    "--restart",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite existing intermediate files.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
def aggregated_classification(input, output, restart, verbose):
    """
    Aggregate the results of the [cyan]marker-classification[/cyan] and
    [cyan]nn-classification[/cyan] modules to classify the sequences in the
    [u]INPUT[/u] file (FASTA format) and write the results to the [u]OUTPUT[/u]
    directory.
    """
    genomad.aggregated_classification.main(input, output, restart, verbose)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.version_option(prog_name="geNomad")
@click.option(
    "--composition",
    type=click.Choice(["auto", "metagenome", "virome"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="Method for estimating sample composition.",
)
@click.option(
    "--force-auto",
    is_flag=True,
    default=False,
    show_default=True,
    help="Force automatic composition estimation regardless of the sample size.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
def score_calibration(input, output, composition, force_auto, verbose):
    """
    Performs score calibration of the sequences in the [u]INPUT[/u] file (FASTA
    format) using the batch correction method and write the results to the
    [u]OUTPUT[/u] directory. This module requires that at least one of the
    classification modules was executed previously ([cyan]marker-classification[/cyan],
    [cyan]nn-classification[/cyan], [cyan]aggregated-classification[/cyan]).
    """
    genomad.score_calibration.main(input, output, composition, force_auto, verbose)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.version_option(prog_name="geNomad")
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
@click.option(
    "--min-score",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0,
    show_default=True,
    help="""Minimum score to flag a sequence as virus or plasmid. By default,
            the sequence is classified as virus/plasmid if its virus/plasmid
            score is higher than its chromosome score, regardless of the value.""",
)
@click.option(
    "--max-fdr",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0.1,
    show_default=True,
    help="""Maximum accepted false discovery rate. This option will be ignored
            if the scores were not calibrated.""",
)
@click.option(
    "--min-plasmid-marker-enrichment",
    type=float,
    default=0.0,
    show_default=True,
    help="""Minimum allowed value for the plasmid marker enrichment score, which
            represents the total enrichment of plasmid markers in the sequence.
            Sequences with multiple plasmid markers will have higher values than
            the ones that encode few or no markers. This option will be ignored
            if the annotation module was not executed.""",
)
@click.option(
    "--min-virus-marker-enrichment",
    type=float,
    default=0.0,
    show_default=True,
    help="""Minimum allowed value for the virus marker enrichment score, which
            represents the total enrichment of virus markers in the sequence.
            Sequences with multiple virus markers will have higher values than
            the ones that encode few or no markers. This option will be ignored
            if the annotation module was not executed.""",
)
@click.option(
    "--min-plasmid-hallmarks",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="""Minimum number of plasmid hallmarks in the identified plasmids. This
            option will be ignored if the annotation module was not executed.""",
)
@click.option(
    "--min-virus-hallmarks",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="""Minimum number of virus hallmarks in the identified viruses. This
            option will be ignored if the annotation module was not executed.""",
)
@click.option(
    "--max-uscg",
    type=int,
    default=4,
    show_default=True,
    help="""Maximum allowed number of universal single copy genes (USCGs) in a
            virus or a plasmid. Sequences with more than this number of USCGs will
            not be classified as viruses or plasmids, regardless of their score.
            This option will be ignored if the annotation module was not executed.""",
)
def summary(
    input,
    output,
    verbose,
    min_score,
    max_fdr,
    min_plasmid_marker_enrichment,
    min_virus_marker_enrichment,
    min_plasmid_hallmarks,
    min_virus_hallmarks,
    max_uscg,
):
    """
    Generates a classification report file for the sequences in the [u]INPUT[/u]
    file (FASTA format) and write it to the [u]OUTPUT[/u] directory. This module
    requires that at least one of the base classification modules was executed
    previously ([cyan]marker-classification[/cyan], [cyan]nn-classification[/cyan]).
    """
    genomad.summary.main(
        input,
        output,
        verbose,
        min_score,
        max_fdr,
        min_plasmid_marker_enrichment,
        min_virus_marker_enrichment,
        min_plasmid_hallmarks,
        min_virus_hallmarks,
        max_uscg,
    )


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.argument("database", type=click.Path(path_type=Path, exists=True))
@click.version_option(prog_name="geNomad")
@click.option(
    "--restart",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite existing intermediate files.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="Number of threads to use.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
@click.option(
    "--cleanup",
    is_flag=True,
    default=False,
    show_default=True,
    help="Delete intermediate files after execution.",
)
@click.option(
    "--disable-find-proviruses",
    is_flag=True,
    default=False,
    show_default=True,
    help="Skip the execution of the [cyan]find-proviruses[/cyan] module.",
)
@click.option(
    "--disable-nn-classification",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Skip the execution of the [cyan]nn-classification[/cyan] and
            [cyan]aggregated-classification[/cyan] modules.""",
)
@click.option(
    "--enable-score-calibration",
    is_flag=True,
    default=False,
    show_default=True,
    help="Execute the [cyan]score-calibration[/cyan] module.",
)
@click.option(
    "--sensitivity",
    "-s",
    type=click.FloatRange(min=0.0),
    default=4.0,
    show_default=True,
    help="""MMseqs2 marker search sensitivity. Higher values will annotate more
            proteins, but the search will be slower and consume more memory.""",
)
@click.option(
    "--splits",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="""Split the data for the MMseqs2 search. Higher values will reduce
            memory usage, but will make the search slower. If the MMseqs2 search
            is failing, try to increase the number of splits.""",
)
@click.option(
    "--skip-integrase-identification",
    is_flag=True,
    default=False,
    show_default=True,
    help="Disable provirus boundary extension using nearby integrases.",
)
@click.option(
    "--skip-trna-identification",
    is_flag=True,
    default=False,
    show_default=True,
    help="Disable provirus boundary extension using nearby tRNAs.",
)
@click.option(
    "--composition",
    type=click.Choice(["auto", "metagenome", "virome"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="Method for estimating sample composition.",
)
@click.option(
    "--force-auto",
    is_flag=True,
    default=False,
    show_default=True,
    help="Force automatic composition estimation regardless of the sample size.",
)
@click.option(
    "--min-score",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0,
    show_default=True,
    help="""Minimum score to flag a sequence as virus or plasmid. By default,
            the sequence is classified as virus/plasmid if its virus/plasmid
            score is higher than its chromosome score, regardless of the value.""",
)
@click.option(
    "--max-fdr",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0.1,
    show_default=True,
    help="""Maximum accepted false discovery rate. This option will be ignored
            if the scores were not calibrated.""",
)
@click.option(
    "--min-plasmid-marker-enrichment",
    type=float,
    default=0.0,
    show_default=True,
    help="""Minimum allowed value for the plasmid marker enrichment score, which
            represents the total enrichment of plasmid markers in the sequence.
            Sequences with multiple plasmid markers will have higher values than
            the ones that encode few or no markers. This option will be ignored
            if the annotation module was not executed.""",
)
@click.option(
    "--min-virus-marker-enrichment",
    type=float,
    default=0.0,
    show_default=True,
    help="""Minimum allowed value for the virus marker enrichment score, which
            represents the total enrichment of virus markers in the sequence.
            Sequences with multiple virus markers will have higher values than
            the ones that encode few or no markers. This option will be ignored
            if the annotation module was not executed.""",
)
@click.option(
    "--min-plasmid-hallmarks",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="""Minimum number of plasmid hallmarks in the identified plasmids. This
            option will be ignored if the annotation module was not executed.""",
)
@click.option(
    "--min-virus-hallmarks",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="""Minimum number of virus hallmarks in the identified viruses. This
            option will be ignored if the annotation module was not executed.""",
)
@click.option(
    "--max-uscg",
    type=int,
    default=4,
    show_default=True,
    help="""Maximum allowed number of universal single copy genes (USCGs) in a
            virus or a plasmid. Sequences with more than this number of USCGs will
            not be classified as viruses or plasmids, regardless of their score.
            This option will be ignored if the annotation module was not executed.""",
)
@click.pass_context
def end_to_end(
    ctx,
    input,
    output,
    database,
    restart,
    threads,
    cleanup,
    verbose,
    disable_find_proviruses,
    disable_nn_classification,
    enable_score_calibration,
    sensitivity,
    splits,
    skip_integrase_identification,
    skip_trna_identification,
    composition,
    force_auto,
    min_score,
    max_fdr,
    min_plasmid_marker_enrichment,
    min_virus_marker_enrichment,
    min_plasmid_hallmarks,
    min_virus_hallmarks,
    max_uscg,
):
    """
    Takes an [u]INPUT[/u] file (FASTA format) and executes all modules of the
    geNomad pipeline for plasmid and virus identification. Output files are
    written in the [u]OUTPUT[/u] directory. A local copy of geNomad's database
    ([u]DATABASE[/u] directory), which can be downloaded with the
    [cyan]download-database[/cyan] command, is required. [yellow]The
    [cyan]end-to-end[/cyan] command omits some options. If you want to have
    a more granular control over the execution parameters, please execute each
    module separately.[/yellow]
    \b
    \n
    \b
    \n
    ╭── geNomad pipeline ────────────────────────────────────╮\n
    │                    ╭───────────╮                       │\n
    │           ╭────────┤   [u]INPUT[/u]   │                       │\n
    │           │        ╰─────┬─────╯                       │\n
    │           │  ╭───────────▼────────────╮                │\n
    │           │  │        annotate        ├──┬──────────╮  │\n
    │           │  ╰───────────┬────────────╯  │          │  │\n
    │           │  ╭───────────▼────────────╮  │          │  │\n
    │           │  │    find-proviruses     │  │          │  │\n
    │           │  ╰┬──────────────────────┬╯  │          │  │\n
    │  ╭────────▼───▼──────────╮╭──────────▼───▼────────╮ │  │\n
    │  │   nn-classification   ││ marker-classification │ │  │\n
    │  ╰────────────┬──────────╯╰──────────┬────────────╯ │  │\n
    │            ╭──▼──────────────────────▼──╮           │  │\n
    │            │ aggregated-classification  │           │  │\n
    │            ╰─────────────┬──────────────╯           │  │\n
    │            ╭ ─ ─ ─ ─ ─ ─ ▼ ─ ─ ─ ─ ─ ─ ─            │  │\n
    │                  score-calibration      │           │  │\n
    │            │ (optional, off by default)             │  │\n
    │             ─ ─ ─ ─ ─ ─ ─┬─ ─ ─ ─ ─ ─ ─ ╯           │  │\n
    │            ╭─────────────▼──────────────╮           │  │\n
    │            │          summary           ◀───────────╯  │\n
    │            ╰─────────────┬──────────────╯              │\n
    │                    ╭─────▼─────╮                       │\n
    │                    │  [u]OUTPUTS[/u]  │                       │\n
    │                    ╰───────────╯                       │\n
    ╰────────────────────────────────────────────────────────╯\n
    """
    ctx.invoke(
        annotate,
        input=input,
        output=output,
        database=database,
        restart=restart,
        threads=threads,
        verbose=verbose,
        cleanup=cleanup,
        sensitivity=sensitivity,
        splits=splits,
    )
    if not disable_find_proviruses:
        ctx.invoke(
            find_proviruses,
            input=input,
            output=output,
            database=database,
            skip_integrase_identification=skip_integrase_identification,
            skip_trna_identification=skip_trna_identification,
            cleanup=cleanup,
            verbose=verbose,
            restart=restart,
            threads=threads,
        )
    ctx.invoke(
        marker_classification,
        input=input,
        output=output,
        database=database,
        restart=restart,
        threads=threads,
        verbose=verbose,
    )
    if not disable_nn_classification:
        ctx.invoke(
            nn_classification,
            input=input,
            output=output,
            cleanup=cleanup,
            verbose=verbose,
            restart=restart,
        )
        ctx.invoke(
            aggregated_classification,
            input=input,
            output=output,
            verbose=verbose,
            restart=restart,
        )
    if enable_score_calibration:
        ctx.invoke(
            score_calibration,
            input=input,
            output=output,
            composition=composition,
            force_auto=force_auto,
            verbose=verbose,
        )
    ctx.invoke(
        summary,
        input=input,
        output=output,
        min_score=min_score,
        max_fdr=max_fdr,
        min_plasmid_marker_enrichment=min_plasmid_marker_enrichment,
        min_virus_marker_enrichment=min_virus_marker_enrichment,
        min_plasmid_hallmarks=min_plasmid_hallmarks,
        min_virus_hallmarks=min_virus_hallmarks,
        max_uscg=max_uscg,
        verbose=verbose,
    )


if __name__ == "__main__":
    cli()
