import math
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from genomad import sequence, utils


class Prodigal:
    def __init__(self, input_file: Path, prodigal_output: Path) -> None:
        self.input_file = input_file
        self.prodigal_output = prodigal_output

    def _run_prodigal_chunk(
        self, workdir_path: Path, current_chunk: int, chunk_file: Path
    ) -> None:
        output_file = workdir_path / f"chunk{current_chunk}.faa"
        cmd = [
            "prodigal-gv",
            "-q",
            "-m",
            "-p",
            "meta",
            "-i",
            chunk_file,
            "-a",
            output_file,
        ]
        subprocess.run(cmd, shell=False, check=True, stdout=subprocess.DEVNULL)

    def _append_prodigal_fasta(self, filepath: Path, start_num: int) -> None:
        pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
        with open(self.prodigal_output, "a") as fout:
            for line in utils.read_file(filepath):
                if line.startswith(">"):
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        start_num += 1
                    line = (
                        match.group(1)
                        + str(start_num)
                        + "_"
                        + match.group(3)
                        + match.group(4)
                    )
                fout.write(line.strip() + "\n")

    def run_parallel_prodigal(self, threads: int) -> None:
        n_sequences = sequence.count_seqs(self.input_file)
        seqs_per_chunk = math.ceil(n_sequences / threads)
        # Delete output file if it already exists
        if self.prodigal_output.is_file():
            self.prodigal_output.unlink()
        seqcnt = 0
        current_chunk = 1
        workdir = tempfile.TemporaryDirectory()
        workdir_path = Path(workdir.name)
        executor = ThreadPoolExecutor(max_workers=threads)
        current_file_path = workdir_path / f"chunk{current_chunk}.fna"
        current_file = open(current_file_path, "w")
        for line in utils.read_file(self.input_file):
            if line[0] == ">" and seqcnt == seqs_per_chunk:
                current_file.close()
                executor.submit(
                    self._run_prodigal_chunk,
                    workdir_path,
                    current_chunk,
                    current_file_path,
                )
                current_file = None
                seqcnt = 0
                current_chunk += 1
            if current_file is None:
                current_file_path = workdir_path / f"chunk{current_chunk}.fna"
                current_file = open(current_file_path, "w")
            current_file.write(line)
            if line[0] == ">":
                seqcnt += 1
        if seqcnt > 0:
            current_file.close()
            executor.submit(
                self._run_prodigal_chunk, workdir_path, current_chunk, current_file_path
            )
        # await completion of tasks
        executor.shutdown(wait=True)
        # collect output
        protid_start = 0
        for c in range(1, current_chunk + 1):
            current_file_path = workdir_path / f"chunk{c}.faa"
            self._append_prodigal_fasta(current_file_path, protid_start)
            protid_start += seqs_per_chunk

    def proteins(self):
        header_parser = re.compile(
            r"(.+)_(.+) # ([0-9]+) # ([0-9]+) # (-1|1) .+rbs_motif=(.+?)"
            r";.+;genetic_code=(.+?);gc_cont=(.+)"
        )
        if not self.prodigal_output.is_file():
            raise FileNotFoundError(f"{self.prodigal_output} was not found.")
        for seq in sequence.read_fasta(self.prodigal_output):
            contig, gene, start, end, strand, rbs, code, gc = header_parser.match(
                seq.header
            ).groups()
            yield (
                contig,
                gene,
                int(start),
                int(end),
                int(strand),
                rbs,
                int(code),
                float(gc),
            )
