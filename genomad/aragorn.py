import math
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from genomad import sequence, utils


class Aragorn:
    def __init__(
        self, input_file: Path, aragorn_output: Path, score_threshold: float = 1.05
    ) -> None:
        self.input_file = input_file
        self.aragorn_output = aragorn_output
        self.score_threshold = int(score_threshold * 100)

    def _run_aragorn_chunk(
        self, workdir_path: Path, current_chunk: int, chunk_file: Path
    ) -> None:
        output_file = workdir_path / f"chunk{current_chunk}.txt"
        cmd = [
            "aragorn",
            "-l",
            f"-ps{self.score_threshold}",
            "-w",
            "-o",
            output_file,
            chunk_file,
        ]
        subprocess.run(cmd, shell=False, check=True, stdout=subprocess.DEVNULL)

    def _append_aragorn_tsv(self, filepath: Path) -> None:
        pattern = re.compile(r"[0-9]+\s+?tRNA-(\S+)\s+?c\[([0-9]+),([0-9]+)\].+")
        current_intervals = []
        with open(self.aragorn_output, "a") as fout, open(filepath) as fin:
            for line in fin:
                if line.startswith(">end"):
                    break
                elif line.startswith(">"):
                    if len(current_intervals):
                        for i, (a, s, e) in enumerate(current_intervals, 1):
                            fout.write(f"{current_contig}_tRNA{i}_{a}\t{s}\t{e}\n")
                    current_contig = line[1:].strip().split()[0]
                    current_intervals = []
                    _ = next(fin)
                elif m := pattern.search(line):
                    a, *coordinates = m.groups()
                    interval = [a, *sorted(map(int, coordinates))]
                    current_intervals.append(interval)
            for i, (a, s, e) in enumerate(current_intervals, 1):
                fout.write(f"{current_contig}_tRNA{i}_{a}\t{s}\t{e}\n")

    def run_parallel_aragorn(self, threads: int) -> None:
        n_sequences = sequence.count_seqs(self.input_file)
        seqs_per_chunk = math.ceil(n_sequences / threads)
        # Delete output file if it already exists
        if self.aragorn_output.is_file():
            self.aragorn_output.unlink()
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
                    self._run_aragorn_chunk,
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
                self._run_aragorn_chunk, workdir_path, current_chunk, current_file_path
            )
        # await completion of tasks
        executor.shutdown(wait=True)
        # collect output
        for c in range(1, current_chunk + 1):
            current_file_path = workdir_path / f"chunk{c}.txt"
            self._append_aragorn_tsv(current_file_path)
