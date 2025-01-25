import multiprocessing.pool
import re
from pathlib import Path

import pyrodigal_gv
from genomad import sequence


_orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True, mask=True)


def _predict_genes(seq):
    return (seq.accession, _orf_finder.find_genes(seq.seq))


class Prodigal:
    def __init__(self, input_file: Path, prodigal_output: Path) -> None:
        self.input_file = input_file
        self.prodigal_output = prodigal_output

    def run_parallel_prodigal(self, threads: int) -> None:
        input_sequences = sequence.read_fasta(self.input_file)
        with (
            multiprocessing.pool.Pool(threads) as pool,
            open(self.prodigal_output, "w") as fout,
        ):
            for seq_i, (seq_acc, predicted_genes) in enumerate(
                pool.imap(_predict_genes, input_sequences), 1
            ):
                for gene_i, gene in enumerate(predicted_genes, 1):
                    header = (
                        f"{seq_acc}_{gene_i} # {gene.begin} # {gene.end} # "
                        + f"{gene.strand} # ID={seq_i}_{gene_i};"
                        + f"partial={int(gene.partial_begin)}{int(gene.partial_end)};"
                        + f"start_type={gene.start_type};rbs_motif={gene.rbs_motif};"
                        + f"rbs_spacer={gene.rbs_spacer};"
                        + f"genetic_code={gene.translation_table};"
                        + f"gc_cont={gene.gc_cont:.3f}"
                    )
                    gene = sequence.Sequence(header, gene.translate(include_stop=False))
                    fout.write(str(gene))

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
