from __future__ import annotations

import gzip
import textwrap
from pathlib import Path
from typing import List

from genomad import utils
from numba import njit


class Sequence:
    def __init__(self, header: str, seq: str, compress: bool = False) -> None:
        self._compress = compress
        self._header = header
        if self._compress:
            self._seq = gzip.compress(seq.encode("ascii"), 1)
        else:
            self._seq = seq.encode("ascii")

    @property
    def header(self) -> str:
        return self._header

    @property
    def accession(self) -> str:
        return self._header.split()[0]

    @property
    def seq(self) -> str:
        if self._compress:
            return gzip.decompress(self._seq).decode()
        else:
            return self._seq.decode()

    @property
    def seq_ascii(self) -> bytes:
        return self.seq.upper().encode("ascii")

    def count(self, substring: str) -> int:
        return self.seq.count(substring)

    def rc(self) -> Sequence:
        tab = self.seq.maketrans("ACTGNactgn", "TGACNtgacn")
        return Sequence(self.header, self.seq.translate(tab)[::-1], self._compress)

    def has_dtr(self, min_length: int = 21) -> bool:
        substring = self.seq.upper()[:min_length]
        pos = self.seq.upper().rfind(substring)
        if pos < len(self) / 2:
            return False
        substring = self.seq.upper()[pos:]
        return self.seq.upper()[: len(substring)] == substring

    def has_itr(self, min_len: int = 21) -> bool:
        rev = self.rc().seq
        return self.seq.upper()[:min_len] == rev.upper()[:min_len]

    def __str__(self) -> str:
        return f">{self.header}\n{textwrap.fill(self.seq, 60)}\n"

    def __repr__(self) -> str:
        if len(self) > 40:
            start = self.seq[:34]
            end = self.seq[-3:]
            seq = f"{start}...{end}"
        else:
            seq = self.seq
        return f"Sequence({self.accession}, {seq})"

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, k: int) -> Sequence:
        return Sequence(self.header, self.seq[k], self._compress)

    def __eq__(self, other: object) -> bool:
        if other.__class__ is self.__class__:
            return self.seq.upper() == other.seq.upper()
        elif other.__class__ is str:
            return self.seq.upper() == other.upper()
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.seq.upper())

    def __add__(self, other: object) -> Sequence:
        if other.__class__ is not self.__class__:
            return NotImplemented
        compress = other._compress or self._compress
        return Sequence(
            f"{self.accession}+{other.accession}", f"{self.seq}{other.seq}", compress
        )


def read_fasta(filepath, uppercase=False, strip_n=False, compress=False):
    with utils.open_file(filepath) as fin:
        last = None
        while True:
            if not last:
                for l in fin:
                    if l[0] == ">":
                        last = l[:-1]
                        break
            if not last:
                break
            name, seqs, last = last[1:], [], None
            for l in fin:
                if l[0] == ">":
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            seqs = "".join(seqs)
            if uppercase:
                seqs = seqs.upper()
            if strip_n:
                seqs = seqs.strip("nN")
            if len(seqs):
                yield Sequence(name, seqs, compress)
            if not last:
                break


def check_fasta(filepath):
    seq_accessions = [seq.accession for seq in read_fasta(filepath)]
    if not len(seq_accessions):
        return False
    elif len(seq_accessions) > len(set(seq_accessions)):
        return False
    else:
        return True


def count_seqs(filepath: Path) -> int:
    return sum(line.startswith(">") for line in utils.read_file(filepath))


def filter_fasta(
    input_filepath, output_filepath, selected_seqs, ignore_gene_suffix=False
):
    with open(output_filepath, "w") as fout:
        for seq in read_fasta(input_filepath):
            seq_name = (
                seq.accession.rsplit("_", 1)[0] if ignore_gene_suffix else seq.accession
            )
            if seq_name in selected_seqs:
                fout.write(f"{seq}\n")


def seq_windows(
    seq: Sequence,
    length: int,
    min_length: int = 0,
    force_first_window: bool = True,
    max_windows: int = None,
) -> Sequence:
    win = 0
    while win * length < len(seq):
        seq_window = seq[win * length : (win + 1) * length]
        if len(seq_window) < min_length:
            if win == 0 and force_first_window:
                yield seq_window
            break
        yield seq_window
        win += 1
        if max_windows and win == max_windows:
            break


@njit
def tokenize_dna(seq: bytes, word_size: int) -> List[int]:
    final_length = len(seq) - word_size + 1
    tokenized_seq = [int(i) for i in range(0)]
    kmer = 0
    countdown = word_size - 1
    mask = (1 << 2 * word_size) - 1
    for base in seq:
        if base == 65:
            kmer = ((kmer << 2) | 0) & mask
        elif base == 67:
            kmer = ((kmer << 2) | 1) & mask
        elif base == 71:
            kmer = ((kmer << 2) | 2) & mask
        elif base == 84:
            kmer = ((kmer << 2) | 3) & mask
        else:
            tokenized_seq += [0] * (word_size - countdown)
            countdown = word_size
        if countdown == 0:
            tokenized_seq.append(kmer + 1)
        else:
            countdown -= 1
    return tokenized_seq[:final_length]
