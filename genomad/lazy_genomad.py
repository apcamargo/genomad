import os
import rich_click as click
from importlib import resources
from genomad.modules.lazy_group import LazyGroup #, help_long 

@click.group(name="Main",
             cls=LazyGroup,
             context_settings={'show_default': True, "help_option_names": ['-h',"-H", '--help']},
             lazy_subcommands={
                 "download-database": "genomad.modules.download.download_database",
                 "end-to-end": "genomad.modules.endtoend.end_to_end",
                 "annotate": "genomad.modules.annotate.annotate",
                 "find-proviruses": "genomad.modules.find_proviruses.find_proviruses",
                 "marker-classification": "genomad.modules.marker_classification.marker_classification",
                 "summary": "genomad.modules.summary.summary",
                 "nn-classification": "genomad.modules.nn_classification.nn_classification",
                #  "aggregated-classification": "genomad.modules.aggregated_classification.aggregated_classification", #  where is this
                #  "plasmid-score": "genomad.modules.plasmid_score.plasmid_score",
                #  "virus-score": "genomad.modules.virus_score.virus_score",
                 "taxonomy": "genomad.modules.taxonomy.taxonomy",
                 "mini-annotate": "genomad.modules.mini_annotate.mini_annotate",
                 "convert-genbank": "genomad.modules.convert_genbank.convert_genbank",
                 "convert-fasta": "genomad.modules.convert_fasta.convert_fasta",
                 "help": "genomad.modules.lazy_group.help_long",
             }
             )
@click.version_option(prog_name="geNomad")
def genomad():
    """geNomad: Identification of mobile genetic elements"""
    pass

if __name__ == "__main__":
    genomad()
