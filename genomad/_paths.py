from dataclasses import dataclass
from pathlib import Path


class classproperty(property):
    def __get__(self, cls, owner):
        return classmethod(self.fget).__get__(None, owner)()


@dataclass
class GenomadData:
    data_dir: Path = Path(__file__).parents[0] / "data"

    @classproperty
    def decision_forest_file(cls) -> Path:
        return cls.data_dir / "decision_forest.ubj"

    @classproperty
    def nn_model_file(cls) -> Path:
        return cls.data_dir / "nn_classifier.h5"

    @classproperty
    def provirus_tagger_file(cls) -> Path:
        return cls.data_dir / "provirus_tagger.crfsuite"

    @classproperty
    def rbs_file(cls) -> Path:
        return cls.data_dir / "rbs_categories.tsv"

    @classproperty
    def score_calibration_weights_file(cls) -> Path:
        return cls.data_dir / "score_calibration_weights.npz"


@dataclass
class GenomadOutputs:
    prefix: str
    output_dir: Path

    # annotate

    @property
    def annotate_log(self) -> Path:
        return self.output_dir / f"{self.prefix}_annotate.log"

    @property
    def annotate_dir(self) -> Path:
        return self.output_dir / f"{self.prefix}_annotate"

    @property
    def annotate_execution_info(self) -> Path:
        return self.annotate_dir / f"{self.prefix}_annotate.json"

    @property
    def annotate_proteins_output(self) -> Path:
        return self.annotate_dir / f"{self.prefix}_proteins.faa"

    @property
    def annotate_mmseqs2_dir(self) -> Path:
        return self.annotate_dir / f"{self.prefix}_mmseqs2"

    @property
    def annotate_mmseqs2_output(self) -> Path:
        return self.annotate_dir / f"{self.prefix}_mmseqs2.tsv"

    @property
    def annotate_genes_output(self) -> Path:
        return self.annotate_dir / f"{self.prefix}_genes.tsv"

    @property
    def annotate_taxonomy_output(self) -> Path:
        return self.annotate_dir / f"{self.prefix}_taxonomy.tsv"

    # find-proviruses

    @property
    def find_proviruses_log(self) -> Path:
        return self.output_dir / f"{self.prefix}_find_proviruses.log"

    @property
    def find_proviruses_dir(self) -> Path:
        return self.output_dir / f"{self.prefix}_find_proviruses"

    @property
    def find_proviruses_execution_info(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_find_proviruses.json"

    @property
    def find_proviruses_output(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus.tsv"

    @property
    def find_proviruses_genes_output(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus_genes.tsv"

    @property
    def find_proviruses_proteins_output(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus_proteins.faa"

    @property
    def find_proviruses_nucleotide_output(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus.fna"

    @property
    def find_proviruses_mmseqs2_input(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus_mmseqs2_input.faa"

    @property
    def find_proviruses_mmseqs2_dir(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus_mmseqs2"

    @property
    def find_proviruses_mmseqs2_output(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus_mmseqs2.tsv"

    @property
    def find_proviruses_aragorn_input(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus_aragorn_input.fna"

    @property
    def find_proviruses_aragorn_output(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus_aragorn.tsv"

    # marker-classification

    @property
    def marker_classification_log(self) -> Path:
        return self.output_dir / f"{self.prefix}_marker_classification.log"

    @property
    def marker_classification_dir(self) -> Path:
        return self.output_dir / f"{self.prefix}_marker_classification"

    @property
    def marker_classification_execution_info(self) -> Path:
        return (
            self.marker_classification_dir / f"{self.prefix}_marker_classification.json"
        )

    @property
    def features_output(self) -> Path:
        return self.marker_classification_dir / f"{self.prefix}_features.tsv"

    @property
    def features_npz_output(self) -> Path:
        return self.marker_classification_dir / f"{self.prefix}_features.npz"

    @property
    def marker_classification_output(self) -> Path:
        return (
            self.marker_classification_dir / f"{self.prefix}_marker_classification.tsv"
        )

    @property
    def marker_classification_npz_output(self) -> Path:
        return (
            self.marker_classification_dir / f"{self.prefix}_marker_classification.npz"
        )

    @property
    def provirus_features_output(self) -> Path:
        return self.marker_classification_dir / f"{self.prefix}_provirus_features.tsv"

    @property
    def provirus_features_npz_output(self) -> Path:
        return self.marker_classification_dir / f"{self.prefix}_provirus_features.npz"

    @property
    def provirus_marker_classification_output(self) -> Path:
        return (
            self.marker_classification_dir
            / f"{self.prefix}_provirus_marker_classification.tsv"
        )

    @property
    def provirus_marker_classification_npz_output(self) -> Path:
        return (
            self.marker_classification_dir
            / f"{self.prefix}_provirus_marker_classification.npz"
        )

    @property
    def find_proviruses_taxonomy_output(self) -> Path:
        return self.find_proviruses_dir / f"{self.prefix}_provirus_taxonomy.tsv"

    # nn-classification

    @property
    def nn_classification_log(self) -> Path:
        return self.output_dir / f"{self.prefix}_nn_classification.log"

    @property
    def nn_classification_dir(self) -> Path:
        return self.output_dir / f"{self.prefix}_nn_classification"

    @property
    def nn_classification_execution_info(self) -> Path:
        return self.nn_classification_dir / f"{self.prefix}_nn_classification.json"

    @property
    def encoded_sequences_dir(self) -> Path:
        return self.nn_classification_dir / f"{self.prefix}_encoded_sequences"

    @property
    def seq_window_id_output(self) -> Path:
        return self.encoded_sequences_dir / f"{self.prefix}_seq_window_id.npz"

    @property
    def nn_classification_output(self) -> Path:
        return self.nn_classification_dir / f"{self.prefix}_nn_classification.tsv"

    @property
    def nn_classification_npz_output(self) -> Path:
        return self.nn_classification_dir / f"{self.prefix}_nn_classification.npz"

    @property
    def encoded_proviruses_dir(self) -> Path:
        return self.nn_classification_dir / f"{self.prefix}_encoded_proviruses"

    @property
    def provirus_window_id_output(self) -> Path:
        return self.encoded_proviruses_dir / f"{self.prefix}_provirus_window_id.npz"

    @property
    def provirus_nn_classification_output(self) -> Path:
        return (
            self.nn_classification_dir / f"{self.prefix}_provirus_nn_classification.tsv"
        )

    @property
    def provirus_nn_classification_npz_output(self) -> Path:
        return (
            self.nn_classification_dir / f"{self.prefix}_provirus_nn_classification.npz"
        )

    # aggregated-classification

    @property
    def aggregated_classification_log(self) -> Path:
        return self.output_dir / f"{self.prefix}_aggregated_classification.log"

    @property
    def aggregated_classification_dir(self) -> Path:
        return self.output_dir / f"{self.prefix}_aggregated_classification"

    @property
    def aggregated_classification_execution_info(self) -> Path:
        return (
            self.aggregated_classification_dir
            / f"{self.prefix}_aggregated_classification.json"
        )

    @property
    def aggregated_classification_output(self) -> Path:
        return (
            self.aggregated_classification_dir
            / f"{self.prefix}_aggregated_classification.tsv"
        )

    @property
    def aggregated_classification_npz_output(self) -> Path:
        return (
            self.aggregated_classification_dir
            / f"{self.prefix}_aggregated_classification.npz"
        )

    @property
    def provirus_aggregated_classification_output(self) -> Path:
        return (
            self.aggregated_classification_dir
            / f"{self.prefix}_provirus_aggregated_classification.tsv"
        )

    @property
    def provirus_aggregated_classification_npz_output(self) -> Path:
        return (
            self.aggregated_classification_dir
            / f"{self.prefix}_provirus_aggregated_classification.npz"
        )

    # score-calibration

    @property
    def score_calibration_log(self) -> Path:
        return self.output_dir / f"{self.prefix}_score_calibration.log"

    @property
    def score_calibration_dir(self) -> Path:
        return self.output_dir / f"{self.prefix}_score_calibration"

    @property
    def score_calibration_execution_info(self) -> Path:
        return self.score_calibration_dir / f"{self.prefix}_score_calibration.json"

    @property
    def score_calibration_compositions_output(self) -> Path:
        return self.score_calibration_dir / f"{self.prefix}_compositions.tsv"

    @property
    def score_calibration_compositions_npz_output(self) -> Path:
        return self.score_calibration_dir / f"{self.prefix}_compositions.npz"

    @property
    def calibrated_marker_classification_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_calibrated_marker_classification.tsv"
        )

    @property
    def calibrated_marker_classification_npz_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_calibrated_marker_classification.npz"
        )

    @property
    def calibrated_nn_classification_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_calibrated_nn_classification.tsv"
        )

    @property
    def calibrated_nn_classification_npz_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_calibrated_nn_classification.npz"
        )

    @property
    def calibrated_aggregated_classification_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_calibrated_aggregated_classification.tsv"
        )

    @property
    def calibrated_aggregated_classification_npz_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_calibrated_aggregated_classification.npz"
        )

    @property
    def provirus_calibrated_marker_classification_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_provirus_calibrated_marker_classification.tsv"
        )

    @property
    def provirus_calibrated_marker_classification_npz_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_provirus_calibrated_marker_classification.npz"
        )

    @property
    def provirus_calibrated_nn_classification_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_provirus_calibrated_nn_classification.tsv"
        )

    @property
    def provirus_calibrated_nn_classification_npz_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_provirus_calibrated_nn_classification.npz"
        )

    @property
    def provirus_calibrated_aggregated_classification_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_provirus_calibrated_aggregated_classification.tsv"
        )

    @property
    def provirus_calibrated_aggregated_classification_npz_output(self) -> Path:
        return (
            self.score_calibration_dir
            / f"{self.prefix}_provirus_calibrated_aggregated_classification.npz"
        )

    # summary

    @property
    def summary_log(self) -> Path:
        return self.output_dir / f"{self.prefix}_summary.log"

    @property
    def summary_dir(self) -> Path:
        return self.output_dir / f"{self.prefix}_summary"

    @property
    def summary_execution_info(self) -> Path:
        return self.summary_dir / f"{self.prefix}_summary.json"

    @property
    def summary_virus_output(self) -> Path:
        return self.summary_dir / f"{self.prefix}_virus_summary.tsv"

    @property
    def summary_virus_sequences_output(self) -> Path:
        return self.summary_dir / f"{self.prefix}_virus.fna"

    @property
    def summary_virus_proteins_output(self) -> Path:
        return self.summary_dir / f"{self.prefix}_virus_proteins.faa"

    @property
    def summary_virus_genes_output(self) -> Path:
        return self.summary_dir / f"{self.prefix}_virus_genes.tsv"

    @property
    def summary_plasmid_output(self) -> Path:
        return self.summary_dir / f"{self.prefix}_plasmid_summary.tsv"

    @property
    def summary_plasmid_sequences_output(self) -> Path:
        return self.summary_dir / f"{self.prefix}_plasmid.fna"

    @property
    def summary_plasmid_proteins_output(self) -> Path:
        return self.summary_dir / f"{self.prefix}_plasmid_proteins.faa"

    @property
    def summary_plasmid_genes_output(self) -> Path:
        return self.summary_dir / f"{self.prefix}_plasmid_genes.tsv"
