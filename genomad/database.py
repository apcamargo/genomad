from pathlib import Path

import taxopy
from genomad import utils


class Database:
    def __init__(self, database_directory: Path) -> None:
        self._directory = database_directory
        with open(database_directory / "version.txt") as fin:
            self._version = float(fin.read().strip())

    @property
    def directory(self) -> Path:
        return self._directory

    @property
    def mmseqs2_db(self) -> Path:
        return self.directory / "genomad_db"

    @property
    def mmseqs2_minimal_db(self) -> Path:
        return self.directory / "genomad_mini_db"

    @property
    def mmseqs2_integrase_db(self) -> Path:
        return self.directory / "genomad_integrase_db"

    @property
    def nodes_dmp(self) -> Path:
        return self.directory / "nodes.dmp"

    @property
    def names_dmp(self) -> Path:
        return self.directory / "names.dmp"

    @property
    def version(self) -> float:
        return self._version

    def get_marker_annotation(self) -> dict:
        """
        Returns a dictionary where the keys are marker names and the values are
        the following:
        (1) Universal single copy gene (USCG)
        (2) Plasmid hallmark
        (3) Virus hallmark
        (4) CONJscan annotation
        (5) AMR annotation
        (6) Functional annotation accessions (Pfam, TIGRFAM, COG, KO)
        (7) Functional description
        """
        marker_annotation = {}
        metadata_file = self.directory / "genomad_marker_metadata.tsv"
        for line in utils.read_file(metadata_file, skip_header=True):
            (
                marker,
                *_,
                uscg,
                plasmid_hallmark,
                virus_hallmark,
                conjscan,
                amr,
                accession,
                description,
                _,
            ) = line.strip().split("\t")
            marker_annotation[marker] = (
                int(uscg != "NA"),
                int(plasmid_hallmark),
                int(virus_hallmark),
                conjscan,
                amr,
                accession,
                description,
            )
        return marker_annotation

    def get_marker_features(self) -> dict:
        """
        Returns a dictionary where the keys are marker names and the values are
        the following:
        (1) Specificity class
        (2) Chromosome SPM
        (3) Plasmid SPM
        (4) Virus SPM
        (5) Giant virus marker
        (6) Universal single copy gene (USCG)
        (7) Plasmid hallmark
        (8) Virus hallmark
        """
        marker_features = {}
        metadata_file = self.directory / "genomad_marker_metadata.tsv"
        for line in utils.read_file(metadata_file, skip_header=True):
            (
                marker,
                _,
                specificity_class,
                _,
                spm_c,
                spm_p,
                spm_v,
                gv_marker,
                *_,
                uscg,
                plasmid_hallmark,
                virus_hallmark,
                _,
                _,
                _,
                _,
                _,
            ) = line.strip().split("\t")
            marker_features[marker] = (
                specificity_class,
                float(spm_c),
                float(spm_p),
                float(spm_v),
                int(gv_marker),
                int(uscg != "NA"),
                int(plasmid_hallmark),
                int(virus_hallmark),
            )
        return marker_features

    def get_taxdb(self) -> taxopy.TaxDb:
        """
        Returns a TaxDb object containing the ICTV taxdump data.
        """
        return taxopy.TaxDb(
            nodes_dmp=self.nodes_dmp, names_dmp=self.names_dmp, keep_files=True
        )
