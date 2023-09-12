import os
import shutil
import subprocess
from pathlib import Path

from genomad import database, utils


class MMseqs2:
    def __init__(
        self,
        mmseqs2_output: Path,
        mmseqs2_directory: Path,
        proteins_output: Path,
        genomad_db: database.Database,
        use_minimal_db: bool = False,
        use_integrase_db: bool = False,
    ) -> None:
        self._mmseqs2_output = mmseqs2_output
        self._mmseqs2_directory = mmseqs2_directory
        self._proteins_output = proteins_output
        self._genomad_db = genomad_db
        if use_integrase_db:
            self._mmseqs2_db = genomad_db.mmseqs2_integrase_db
            self._include_taxid = False
        elif use_minimal_db:
            self._mmseqs2_db = genomad_db.mmseqs2_minimal_db
            self._include_taxid = True
        else:
            self._mmseqs2_db = genomad_db.mmseqs2_db
            self._include_taxid = True

    @property
    def mmseqs2_output(self) -> Path:
        return self._mmseqs2_output

    @property
    def mmseqs2_directory(self) -> Path:
        return self._mmseqs2_directory

    @property
    def proteins_output(self) -> Path:
        return self._proteins_output

    @property
    def mmseqs2_db(self) -> Path:
        return self._mmseqs2_db

    @property
    def include_taxid(self) -> bool:
        return self._include_taxid

    def run_mmseqs2(
        self, threads: int, sensitivity: float, evalue: float, splits: int
    ) -> None:
        if self.mmseqs2_directory.exists():
            shutil.rmtree(self.mmseqs2_directory)
        # Create the MMseqs2 output directory and its subdirectories
        self.mmseqs2_directory.mkdir()
        query_db_dir = self.mmseqs2_directory / "query_db"
        query_db_dir.mkdir()
        search_db_dir = self.mmseqs2_directory / "search_db"
        search_db_dir.mkdir()
        besthit_db_dir = self.mmseqs2_directory / "besthit_db"
        besthit_db_dir.mkdir()
        # Define the query, search, and besthit databases
        query_db = query_db_dir / "query_db"
        prefilter_db = search_db_dir / "prefilter_db"
        swapresults_1_db = search_db_dir / "swapresults_1_db"
        align_1_db = search_db_dir / "align_1_db"
        align_2_db = search_db_dir / "align_2_db"
        swapresults_2_db = search_db_dir / "swapresults_2_db"
        besthit_db = besthit_db_dir / "besthit_db"
        # Define the MMseqs2 commands
        createdb_command = ["mmseqs", "createdb", self.proteins_output, query_db]
        prefilter_command = [
            "mmseqs",
            "prefilter",
            query_db,
            self.mmseqs2_db,
            prefilter_db,
            "--threads",
            str(threads),
            "-s",
            str(sensitivity),
            "--split",
            str(splits),
            "--split-mode",
            "0",
            "--max-seqs",
            "10000000",
            "--min-ungapped-score",
            "25",
            "-k",
            "5",
        ]
        swapresults_1_command = [
            "mmseqs",
            "swapresults",
            query_db,
            self.mmseqs2_db,
            prefilter_db,
            swapresults_1_db,
            "--threads",
            str(threads),
        ]
        align_1_command = [
            "mmseqs",
            "align",
            self.mmseqs2_db,
            query_db,
            swapresults_1_db,
            align_1_db,
            "--threads",
            str(threads),
            "--alignment-mode",
            "1",
            "-e",
            str(evalue),
            "--max-rejected",
            "280",
        ]
        align_2_command = [
            "mmseqs",
            "align",
            self.mmseqs2_db,
            query_db,
            align_1_db,
            align_2_db,
            "--threads",
            str(threads),
            "--alignment-mode",
            "2",
            "-e",
            str(evalue),
            "--cov-mode",
            "2",
            "-c",
            "0.2",
        ]
        swapresults_2_command = [
            "mmseqs",
            "swapresults",
            self.mmseqs2_db,
            query_db,
            align_2_db,
            swapresults_2_db,
            "--threads",
            str(threads),
        ]
        besthit_command = [
            "mmseqs",
            "filterdb",
            swapresults_2_db,
            besthit_db,
            "--extract-lines",
            "1",
        ]
        if self.include_taxid:
            output_columns = "qheader,target,evalue,bits,taxid"
        else:
            output_columns = "qheader,target,evalue,bits"
        convertalis_command = [
            "mmseqs",
            "convertalis",
            query_db,
            self.mmseqs2_db,
            besthit_db,
            self.mmseqs2_output,
            "--format-output",
            output_columns,
            "--threads",
            str(threads),
        ]
        log_file = self.mmseqs2_directory / "mmseqs2.log"
        with open(log_file, "w") as fout:
            # Check if the protein FASTA file is not empty
            if os.stat(self._proteins_output).st_size > 0:
                for command in [
                    createdb_command,
                    prefilter_command,
                    swapresults_1_command,
                    align_1_command,
                    align_2_command,
                    swapresults_2_command,
                    besthit_command,
                    convertalis_command,
                ]:
                    try:
                        subprocess.run(command, stdout=fout, stderr=fout, check=True)
                    except subprocess.CalledProcessError as e:
                        command_str = " ".join([str(i) for i in command])
                        raise Exception(f"'{command_str}' failed.") from e
            else:
                fout.write("No queries")
                open(self.mmseqs2_output, "w").close()

    def get_matches(self) -> dict:
        if not self.mmseqs2_output.is_file():
            raise FileNotFoundError(f"{self.mmseqs2_output} was not found.")
        gene_matches = {}
        for line in utils.read_file(self.mmseqs2_output):
            if self.include_taxid:
                gene, match, evalue, bitscore, taxid = line.strip().split("\t")
                gene = gene.split()[0]
                taxid = "1" if taxid == "0" else taxid
                gene_matches[gene] = (match, float(evalue), int(bitscore), int(taxid))
            else:
                gene, match, evalue, bitscore = line.strip().split("\t")
                gene = gene.split()[0]
                gene_matches[gene] = (match, float(evalue), int(bitscore), 1)
        return gene_matches
