import shutil
import sys
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import numpy as np
import pycrfsuite
from genomad import aragorn, database, mmseqs2, sequence, taxonomy, utils
from genomad._paths import GenomadData, GenomadOutputs


@dataclass
class GeneTable:
    seq_name: str
    starts: List[int] = field(default_factory=list)
    ends: List[int] = field(default_factory=list)
    spm_c: List[float] = field(default_factory=list)
    spm_v: List[float] = field(default_factory=list)
    v_vs_c_score: List[float] = field(default_factory=list)
    c_markers: List[bool] = field(default_factory=list)
    v_markers: List[bool] = field(default_factory=list)
    integrases: List[bool] = field(default_factory=list)
    trna_starts: List[int] = field(default_factory=list)
    trna_ends: List[int] = field(default_factory=list)

    @property
    def n_genes(self) -> int:
        return len(self.starts)

    @property
    def n_c_markers(self) -> int:
        return sum(self.c_markers)

    @property
    def n_v_markers(self) -> int:
        return sum(self.v_markers)

    @property
    def integrase_starts(self) -> int:
        return [s for s, i in zip(self.starts, self.integrases) if i]

    @property
    def integrase_ends(self) -> int:
        return [e for e, i in zip(self.ends, self.integrases) if i]


class ProvirusCRF:
    def __init__(self):
        self.tagger = pycrfsuite.Tagger()
        provirus_tagger_file = GenomadData.provirus_tagger_file
        self.tagger.open(str(provirus_tagger_file))

    def score_provirus_genes(
        self,
        spm_v_array,
        spm_c_array,
    ):
        crf_input = [{"spm_v": i, "spm_c": j} for i, j in zip(spm_v_array, spm_c_array)]
        self.tagger.set(crf_input)
        scores = np.array([self.tagger.marginal("V", i) for i in range(len(crf_input))])
        self.tagger.set([{} for _ in range(len(spm_v_array))])
        background_scores = np.array(
            [self.tagger.marginal("V", i) for i in range(len(crf_input))]
        )
        scores = utils.logistic((scores - background_scores), temperature=0.2)
        return scores


@dataclass
class Provirus:
    seq_name: str
    start: int
    end: int
    n_genes: int
    v_vs_c_score: float
    has_integrase: bool
    integrase_indices: List[int]
    is_edge: bool

    def __len__(self):
        return self.end - self.start + 1

    @property
    def provirus_name(self) -> str:
        return f"{self.seq_name}|provirus_{self.start}_{self.end}"


def yield_gene_tables(
    genes_output: Path,
    database_obj: database.Database,
    integrase_mmseqs_output: Optional[Path] = None,
    aragorn_output: Optional[Path] = None,
):
    marker_features_dict = database_obj.get_marker_features()
    integrase_list = []
    # Create a set containing the gene names of detected integrases
    if integrase_mmseqs_output:
        integrase_list.extend(
            line.strip().split("\t")[0]
            for line in utils.read_file(integrase_mmseqs_output)
        )
    integrase_list = set(integrase_list)
    # Create a dict containing the starts and ends of all the tRNAs of a given contig
    trna_dict = defaultdict(lambda: [[], []])
    if aragorn_output:
        for line in utils.read_file(aragorn_output):
            contig, start, end = line.strip().split("\t")
            contig, start, end = contig.rsplit("_", 2)[0], int(start), int(end)
            trna_dict[contig][0].append(start)
            trna_dict[contig][1].append(end)
    # Read the genes output
    current_contig = None
    current_genetable = None
    for line in utils.read_file(genes_output, skip_header=True):
        gene, start, end, _, _, _, _, _, match, *_ = line.strip().split("\t")
        specificity_class, spm_c, _, spm_v, *_ = marker_features_dict.get(
            match, (None, 0.0, 0.0, 0.0, 0)
        )
        contig, start, end = (
            gene.rsplit("_", 1)[0],
            int(start),
            int(end),
        )
        if contig != current_contig:
            if current_contig:
                yield current_genetable
            current_genetable = GeneTable(contig)
            current_contig = contig
        current_genetable.starts.append(start)
        current_genetable.ends.append(end)
        current_genetable.spm_c.append(spm_c)
        current_genetable.spm_v.append(spm_v)
        current_genetable.v_vs_c_score.append(np.exp(spm_v) - np.exp(spm_c))
        if specificity_class:
            current_genetable.c_markers.append(specificity_class.startswith("C"))
            current_genetable.v_markers.append(specificity_class.startswith("V"))
        else:
            current_genetable.c_markers.append(False)
            current_genetable.v_markers.append(False)
        current_genetable.integrases.append(
            integrase_mmseqs_output and gene in integrase_list
        )
        if aragorn_output:
            current_genetable.trna_starts = trna_dict[current_contig][0]
            current_genetable.trna_ends = trna_dict[current_contig][1]
    if current_genetable:
        yield current_genetable


def tag_provirus_genes(
    provirus_array,
    threshold,
    genetable,
    min_markers_host_island=2,
    min_markers_host_edge=1,
    min_genes_host_island=6,
    min_genes_host_edge=4,
    min_markers_phage_island=1,
    min_markers_phage_edge=1,
    min_genes_phage_island=5,
    min_genes_phage_edge=3,
):
    provirus_array = (provirus_array >= threshold).astype(int)
    # Convert small host regions to phage
    total_count = 0
    count_array, value_array = utils.rle_encode(provirus_array)
    for i in range(len(count_array)):
        if value_array[i] == 0:
            spm_c_array = np.array(
                genetable.spm_c[total_count : total_count + count_array[i]]
            )
            spm_v_array = np.array(
                genetable.spm_v[total_count : total_count + count_array[i]]
            )
            n_c_markers = (spm_c_array > spm_v_array).sum()
            n_v_markers = (spm_v_array > spm_c_array).sum()
            # If not in edge
            if (i < len(count_array) - 1) and (i > 0):
                if (
                    (count_array[i] < min_genes_host_island)
                    or (n_c_markers < min_markers_host_island)
                    or (n_c_markers <= n_v_markers)
                ):
                    value_array[i] = 1
            # If in edge
            elif (
                (count_array[i] < min_genes_host_edge)
                or (n_c_markers < min_markers_host_edge)
                or (n_c_markers <= n_v_markers)
            ):
                value_array[i] = 1
        total_count += count_array[i]
    # Convert small phage regions to host
    total_count = 0
    count_array, value_array = utils.rle_encode(
        utils.rle_decode(count_array, value_array)
    )
    for i in range(len(count_array)):
        if value_array[i] == 1:
            spm_c_array = np.array(
                genetable.spm_c[total_count : total_count + count_array[i]]
            )
            spm_v_array = np.array(
                genetable.spm_v[total_count : total_count + count_array[i]]
            )
            n_c_markers = (spm_c_array > spm_v_array).sum()
            n_v_markers = (spm_v_array > spm_c_array).sum()
            # If not in edge
            if (i < len(count_array) - 1) and (i > 0):
                if (
                    (count_array[i] < min_genes_phage_island)
                    or (n_v_markers < min_markers_phage_island)
                    or (n_v_markers <= n_c_markers)
                ):
                    value_array[i] = 0
            # If in edge
            elif (
                (count_array[i] < min_genes_phage_edge)
                or (n_v_markers < min_markers_phage_edge)
                or (n_v_markers <= n_c_markers)
            ):
                value_array[i] = 0
        total_count += count_array[i]
    return utils.rle_decode(count_array, value_array)


def extend_provirus_edges(provirus_labels, genetable, feature_type, max_dist):
    if feature_type == "integrase":
        feature_starts = genetable.integrase_starts
        feature_ends = genetable.integrase_ends
    elif feature_type == "trna":
        feature_starts = genetable.trna_starts
        feature_ends = genetable.trna_ends
    # If the feature type doesn't match "integrase" ot "trna", return the original labels
    else:
        return provirus_labels
    # If the scaffold has only virus or host genes, return the original labels
    if len(set(provirus_labels)) <= 1:
        return provirus_labels
    # Create a provirus coordinate matrix
    total_count = 0
    count_array, value_array = utils.rle_encode(provirus_labels)
    provirus_coordinates = []
    for count, value in zip(count_array, value_array):
        if value == 1:
            provirus_coordinates.append(
                [genetable.starts[total_count], genetable.ends[total_count + count - 1]]
            )
        total_count += count
    provirus_coordinates = np.array(provirus_coordinates)
    # Create the `n_modifications` variable that will keep track whether a provirus
    # coordinate was modified
    n_modifications = 0
    # If the scaffold does not have at least one feature (integrase or tRNA) and a provirus:
    if not (len(feature_starts) and len(provirus_coordinates)):
        return provirus_labels
    # If the scaffold has at least one feature (integrase or tRNA) and a provirus:
    else:
        # Create a matrix to store feature × provirus distances
        feature_provirus_distances_matrix = []
        # Iterate through features
        for (feature_start, feature_end) in zip(feature_starts, feature_ends):
            feature_provirus_distances = []
            # Iterate through proviruses
            for (provirus_start, provirus_end) in provirus_coordinates:
                # If feature after provirus:
                if feature_start > provirus_end:
                    distance = feature_end - provirus_end
                # If feature before provirus:
                elif feature_end < provirus_start:
                    distance = feature_start - provirus_start
                # If feature is inside the provirus:
                else:
                    distance = 0
                feature_provirus_distances.append(distance)
            feature_provirus_distances_matrix.append(feature_provirus_distances)
        feature_provirus_distances_matrix = np.array(feature_provirus_distances_matrix)
        # For each feature, retrieve the closes provirus
        closest_provirus = np.abs(feature_provirus_distances_matrix).argmin(1)
        # Get the relative distances (positive if the feature if after the provirus and
        # negative if the feature is before the provirus)
        relative_distances = feature_provirus_distances_matrix[
            np.arange(feature_provirus_distances_matrix.shape[0]), closest_provirus
        ]
        # Get the absolute distances
        absolute_distances = np.abs(relative_distances)
        # For each provirus, get the closest feature
        closest_feature = np.abs(feature_provirus_distances_matrix).argmin(0)
        # For each feature, identify if the match is reciprocal (that is, if it is also the closest
        # feature to its closest provirus)
        reciprocal = closest_feature[closest_provirus] == np.arange(
            feature_provirus_distances_matrix.shape[0]
        )
        # Filter out features that are more than `max_dist` from the closest provirus
        closest_provirus = closest_provirus[absolute_distances <= max_dist]
        relative_distances = relative_distances[absolute_distances <= max_dist]
        reciprocal = reciprocal[absolute_distances <= max_dist]
        # If there is at least one feature that is less than `max_dist` to a provirus:
        if len(relative_distances):
            for p, d, r in zip(closest_provirus, relative_distances, reciprocal):
                # Get the coordinates of chromosome markers within the scaffold
                c_marker_starts = np.array(genetable.starts)[genetable.c_markers]
                c_marker_ends = np.array(genetable.ends)[genetable.c_markers]
                # If feature is after the provirus (d > 0), is the closest feature to the provirus,
                # and there are no chromosome markers between the feature and the provirus, extend
                # the provirus upstream and updte `n_modifications`
                if (
                    (d > 0)
                    and r
                    and not np.logical_and(
                        c_marker_starts >= provirus_coordinates[p][1],
                        c_marker_ends <= provirus_coordinates[p][1] + d,
                    ).any()
                ):
                    provirus_coordinates[p][1] += d
                    n_modifications += 1
                # If feature is before the provirus (d < 0), is the closest feature to the provirus,
                # and there are no chromosome markers between the feature and the provirus, extend
                # the provirus downstream and updte `n_modifications`
                elif (
                    (d < 0)
                    and r
                    and not np.logical_and(
                        c_marker_ends <= provirus_coordinates[p][0],
                        c_marker_starts >= provirus_coordinates[p][0] - d,
                    ).any()
                ):
                    provirus_coordinates[p][0] += d
                    n_modifications += 1
        # If none of the provirus coordinates were modified, return the original labels:
        if not n_modifications:
            return provirus_labels
        # If the coordinates of at least one provirus was modified, update the provirus labels:
        else:
            updated_labels = []
            # Iterate through the genes and check if they are contained within a provirus
            for gene_start, gene_end in zip(genetable.starts, genetable.ends):
                in_provirus = [
                    gene_start >= p[0] and gene_end <= p[1]
                    for p in provirus_coordinates
                ]
                updated_labels.append(int(any(in_provirus)))
            return updated_labels


def yield_proviruses(
    genetable: GeneTable,
    provirus_labels: List[int],
    threshold: float,
    in_edge_threshold: float,
    has_integrase_threshold: float,
) -> Provirus:
    total_count = 0
    count_array, value_array = utils.rle_encode(provirus_labels)
    n_islands = len(count_array)
    for i, (count, value) in enumerate(zip(count_array, value_array)):
        if value == 1:
            v_vs_c_score = genetable.v_vs_c_score[total_count : total_count + count]
            v_vs_c_score = sum(v_vs_c_score)
            has_integrase = any(genetable.integrases[total_count : total_count + count])
            in_edge = i in [0, n_islands - 1]
            if (
                (in_edge and v_vs_c_score >= in_edge_threshold)
                or (has_integrase and v_vs_c_score >= has_integrase_threshold)
                or (not in_edge and not has_integrase and v_vs_c_score >= threshold)
            ):
                if has_integrase:
                    integrase_indices = np.arange(len(provirus_labels))[
                        total_count : total_count + count
                    ]
                    integrase_indices = integrase_indices[
                        genetable.integrases[total_count : total_count + count]
                    ]
                    integrase_indices = integrase_indices.tolist()
                else:
                    integrase_indices = []
                yield Provirus(
                    genetable.seq_name,
                    genetable.starts[total_count],
                    genetable.ends[total_count + count - 1],
                    count,
                    v_vs_c_score,
                    has_integrase,
                    integrase_indices,
                    in_edge,
                )
        total_count += count


def main(
    input_path,
    output_path,
    database_path,
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
    # Create `output_path` if it does not exist
    if not output_path.is_dir():
        output_path.mkdir()
    # Define the prefix and the output files
    prefix = input_path.stem
    if utils.is_compressed(input_path) != utils.Compression.uncompressed:
        prefix = prefix.rsplit(".", 1)[0]
    outputs = GenomadOutputs(prefix, output_path)
    # Create the console
    console = utils.HybridConsole(
        output_file=outputs.find_proviruses_log, verbose=verbose
    )
    # Create a dictionary containing the parameters
    parameter_dict = {
        "skip_integrase_identification": skip_integrase_identification,
        "skip_trna_identification": skip_trna_identification,
        "crf_threshold": crf_threshold,
        "marker_threshold": marker_threshold,
        "marker_threshold_integrase": marker_threshold_integrase,
        "marker_threshold_edge": marker_threshold_edge,
        "max_integrase_distance": max_integrase_distance,
        "max_trna_distance": max_trna_distance,
        "sensitivity": sensitivity,
        "evalue": evalue,
    }

    # Display the module header
    output_files_path_list = [
        outputs.find_proviruses_execution_info,
        outputs.find_proviruses_output,
        outputs.find_proviruses_nucleotide_output,
        outputs.find_proviruses_proteins_output,
        outputs.find_proviruses_genes_output,
        outputs.find_proviruses_taxonomy_output,
    ]
    output_files_description_list = [
        "execution parameters",
        "provirus data",
        "provirus nucleotide sequences",
        "provirus protein sequences",
        "provirus gene annotation data",
        "provirus taxonomic assignment",
    ]
    if not skip_integrase_identification:
        output_files_path_list.append(outputs.find_proviruses_mmseqs2_output)
        output_files_description_list.append("MMseqs2 output file")
    if not skip_trna_identification:
        output_files_path_list.append(outputs.find_proviruses_aragorn_output)
        output_files_description_list.append("Aragorn output file")
    utils.display_header(
        console,
        "find-proviruses",
        ("This will find putative proviral regions within the input sequences."),
        outputs.find_proviruses_dir,
        output_files_path_list,
        output_files_description_list,
    )

    # Check if the annotation module was previously executed
    if (
        not outputs.annotate_genes_output.exists()
        or not outputs.annotate_proteins_output.exists()
    ):
        console.error(
            f"[u]{outputs.annotate_genes_output.name}[/u] and [u]{outputs.annotate_proteins_output.name}[/u] "
            "were not found in the output directory. Please execute the [cyan]annotate[/cyan] module "
            "to generate them."
        )
        sys.exit(1)
    # Check if the input FASTA for `annotate` and `find-proviruses` is the same
    elif not utils.compare_executions(
        input_path, {}, outputs.annotate_execution_info, only_md5=True
    ):
        console.error(
            "The input FASTA file is different from the one used in the [cyan]annotate[/cyan] "
            "module. Please execute both modules using the same input"
        )
        sys.exit(1)

    # Check if the required binaries are in the PATH
    required_executables = []
    if not skip_integrase_identification:
        required_executables.append("mmseqs")
    if not skip_trna_identification:
        required_executables.append("aragorn")
    if missing_executables := utils.check_executables(required_executables):
        console.error(
            "These dependencies are missing and need to be installed: [u]"
            + "[/u], [u]".join(missing_executables)
            + "[/u]. The current execution will be terminated."
        )
        sys.exit(1)

    # Check the input FASTA file
    if not sequence.check_fasta(input_path):
        console.error(
            f"[u]{input_path}[/u] is either empty or contains multiple entries "
            "with the same identifier. Please check your input FASTA file and "
            "execute [cyan]genomad find-proviruses[/cyan] again."
        )
        sys.exit(1)

    # Print initial log
    console.log("Executing [cyan]genomad find-proviruses[/cyan].")

    # Check if steps can be skipped
    skip = False
    if (
        outputs.find_proviruses_execution_info.exists()
        and any(path.exists() for path in output_files_path_list)
        and not restart
    ):
        # Check if the previous execution used the same input file and parameters
        if utils.compare_executions(
            input_path,
            parameter_dict,
            outputs.find_proviruses_execution_info,
        ):
            skip = True
            console.log(
                "Previous execution detected. Steps will be skipped "
                "unless their outputs are not found. Use the [cyan]"
                "--restart[/cyan] option to force the execution of all "
                "the steps again."
            )
        else:
            console.log(
                "The input file or the parameters changed since the last "
                "execution. Previous outputs will be overwritten."
            )

    # If necessary, create the subdirectory
    if not outputs.find_proviruses_dir.is_dir():
        console.log(
            f"Creating the [green]{outputs.find_proviruses_dir}[/green] directory."
        )
        outputs.find_proviruses_dir.mkdir()

    # Write the execution data to `execution_info_output`
    utils.write_execution_info(
        "find-proviruses",
        input_path,
        parameter_dict,
        outputs.find_proviruses_execution_info,
    )

    # Initialize the database object
    database_obj = database.Database(database_path)

    # Identify the scaffolds thain contain at least one virus and one chromosome markers
    target_contigs = [
        genetable.seq_name
        for genetable in yield_gene_tables(
            outputs.annotate_genes_output,
            database_obj,
            None,
            None,
        )
        if genetable.n_c_markers and genetable.n_v_markers
    ]
    target_contigs = set(target_contigs)

    # If there are no contigs with at least one virus and one chromosome marker,
    # write empty outputs and end the execution.
    if not len(target_contigs):
        console.log("No potential provirus-carrying sequences were identified.")
        for f in output_files_path_list:
            # Don't overwite the execution ifo file
            if f != outputs.find_proviruses_execution_info:
                open(f, "w").close()
        with open(outputs.find_proviruses_output, "w") as fout:
            fout.write(
                "seq_name\tsource_seq\tstart\tend\tlength\tn_genes\t"
                "v_vs_c_score\tin_seq_edge\tintegrases\n"
            )
        with open(outputs.find_proviruses_genes_output, "w") as fout:
            fout.write(
                "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\t"
                "marker\tevalue\tbitscore\tuscg\tannotation_accessions\tannotation_description\n"
            )
        if cleanup and outputs.find_proviruses_mmseqs2_input.exists():
            outputs.find_proviruses_mmseqs2_input.unlink()
        if cleanup and outputs.find_proviruses_aragorn_input.exists():
            outputs.find_proviruses_aragorn_input.unlink()
        console.log("geNomad find-proviruses finished!", style="yellow")
        return

    # Integrase identification
    if skip and outputs.find_proviruses_mmseqs2_output.exists():
        console.log(
            f"[green]{outputs.find_proviruses_mmseqs2_output.name}[/green] was found. "
            "Skipping integrase search."
        )
    elif not skip_integrase_identification:
        with console.status(
            "Identifying integrases with MMseqs2 and geNomad database "
            f"(v{database_obj.version})."
        ):
            sequence.filter_fasta(
                outputs.annotate_proteins_output,
                outputs.find_proviruses_mmseqs2_input,
                target_contigs,
                ignore_gene_suffix=True,
            )
            mmseqs2_obj = mmseqs2.MMseqs2(
                outputs.find_proviruses_mmseqs2_output,
                outputs.find_proviruses_mmseqs2_dir,
                outputs.find_proviruses_mmseqs2_input,
                database_obj,
                use_integrase_db=True,
            )
            mmseqs2_obj.run_mmseqs2(threads, sensitivity, evalue, 0)
            console.log(
                "Integrases identified with MMseqs2 and geNomad database "
                f"(v{database_obj.version}) were written to "
                f"[green]{outputs.find_proviruses_mmseqs2_output.name}[/green]."
            )
    if cleanup and outputs.find_proviruses_mmseqs2_dir.is_dir():
        console.log(
            f"Deleting [green]{outputs.find_proviruses_mmseqs2_dir.name}[/green]."
        )
        shutil.rmtree(outputs.find_proviruses_mmseqs2_dir)
    if cleanup and outputs.find_proviruses_mmseqs2_input.exists():
        console.log(
            f"Deleting [green]{outputs.find_proviruses_mmseqs2_input.name}[/green]."
        )
        outputs.find_proviruses_mmseqs2_input.unlink()

    # tRNA identification
    if skip and outputs.find_proviruses_aragorn_output.exists():
        console.log(
            f"[green]{outputs.find_proviruses_aragorn_output.name}[/green] was found. "
            "Skipping tRNA identification."
        )
    elif not skip_trna_identification:
        with console.status("Identifying tRNAs with Aragorn."):
            sequence.filter_fasta(
                input_path,
                outputs.find_proviruses_aragorn_input,
                target_contigs,
            )
            aragorn_obj = aragorn.Aragorn(
                outputs.find_proviruses_aragorn_input,
                outputs.find_proviruses_aragorn_output,
            )
            aragorn_obj.run_parallel_aragorn(threads)
            console.log(
                "tRNAs identified with Aragorn were written to "
                f"[green]{outputs.find_proviruses_aragorn_output.name}[/green]."
            )
    if cleanup and outputs.find_proviruses_aragorn_input.exists():
        console.log(
            f"Deleting [green]{outputs.find_proviruses_aragorn_input.name}[/green]."
        )
        outputs.find_proviruses_aragorn_input.unlink()

    # Flag proviral genes using the CRF model
    with console.status("Identifying provirus regions."):
        crf = ProvirusCRF()
        provirus_dict = OrderedDict()
        for genetable in yield_gene_tables(
            outputs.annotate_genes_output,
            database_obj,
            None
            if skip_integrase_identification
            else outputs.find_proviruses_mmseqs2_output,
            None
            if skip_trna_identification
            else outputs.find_proviruses_aragorn_output,
        ):
            if genetable.seq_name in target_contigs:
                scores = crf.score_provirus_genes(genetable.spm_v, genetable.spm_c)
                provirus_labels = tag_provirus_genes(scores, crf_threshold, genetable)
                # Extend boundaries using integrases
                if not skip_integrase_identification:
                    provirus_labels = extend_provirus_edges(
                        provirus_labels, genetable, "integrase", max_integrase_distance
                    )
                # Extend boundaries using tRNAs
                if not skip_trna_identification:
                    provirus_labels = extend_provirus_edges(
                        provirus_labels, genetable, "trna", max_trna_distance
                    )
                # If the label array contains 1s and 0s (both proviral and host regions):
                if len(set(provirus_labels)) > 1:
                    provirus_dict[genetable.seq_name] = list(
                        yield_proviruses(
                            genetable,
                            provirus_labels,
                            threshold=marker_threshold,
                            in_edge_threshold=marker_threshold_edge,
                            has_integrase_threshold=marker_threshold_integrase,
                        )
                    )
        console.log("Provirus regions identified.")

    # Writting outputs
    with console.status(
        f"Writing [green]{outputs.find_proviruses_output.name}[/green]."
    ):
        with open(outputs.find_proviruses_output, "w") as fout:
            fout.write(
                "seq_name\tsource_seq\tstart\tend\tlength\tn_genes\t"
                "v_vs_c_score\tin_seq_edge\tintegrases\n"
            )
            for contig_proviruses in provirus_dict.values():
                for provirus in contig_proviruses:
                    if provirus.has_integrase:
                        integrase_genes = [
                            f"{provirus.provirus_name}_{i + 1}"
                            for i in provirus.integrase_indices
                        ]
                        integrase_genes = ";".join(integrase_genes)
                    else:
                        integrase_genes = "NA"
                    fout.write(
                        f"{provirus.provirus_name}\t"
                        f"{provirus.seq_name}\t"
                        f"{provirus.start}\t"
                        f"{provirus.end}\t"
                        f"{provirus.end - provirus.start + 1}\t"
                        f"{provirus.n_genes}\t"
                        f"{provirus.v_vs_c_score:.4f}\t"
                        f"{provirus.is_edge}\t"
                        f"{integrase_genes}\n"
                    )
        console.log(
            f"Provirus data was written to [green]{outputs.find_proviruses_output.name}[/green]."
        )

    with console.status(
        f"Writing [green]{outputs.find_proviruses_nucleotide_output.name}[/green]."
    ):
        with open(outputs.find_proviruses_nucleotide_output, "w") as fout:
            for seq in sequence.read_fasta(input_path):
                if seq.accession in provirus_dict:
                    for provirus in provirus_dict[seq.accession]:
                        provirus_seq = seq.seq[provirus.start - 1 : provirus.end]
                        provirus_seq = sequence.Sequence(
                            provirus.provirus_name, provirus_seq
                        )
                        fout.write(str(provirus_seq))
        console.log(
            "Provirus nucleotide sequences were written to "
            f"[green]{outputs.find_proviruses_nucleotide_output.name}[/green]."
        )

    with console.status(
        f"Writing [green]{outputs.find_proviruses_proteins_output.name}[/green]."
    ):
        with open(outputs.find_proviruses_proteins_output, "w") as fout:
            for seq in sequence.read_fasta(outputs.annotate_proteins_output):
                contig = seq.accession.rsplit("_", 1)[0]
                if contig in provirus_dict:
                    start = int(seq.header.split()[2])
                    end = int(seq.header.split()[4])
                    in_provirus = np.array(
                        [
                            start >= i.start and end <= i.end
                            for i in provirus_dict[contig]
                        ]
                    )
                    if any(in_provirus):
                        provirus_name = provirus_dict[contig][
                            np.argwhere(in_provirus)[0][0]
                        ].provirus_name
                        gene_number = seq.accession.rsplit("_", 1)[1]
                        gene_name = f"{provirus_name}_{gene_number}"
                        header = f"{gene_name} {seq.header.split(maxsplit=1)[1]}"
                        seq = sequence.Sequence(header, seq.seq)
                        fout.write(str(seq))
        console.log(
            "Provirus protein sequences were written to "
            f"[green]{outputs.find_proviruses_proteins_output.name}[/green]."
        )

    with console.status(
        f"Writing [green]{outputs.find_proviruses_genes_output.name}[/green]."
    ):
        with open(outputs.find_proviruses_genes_output, "w") as fout:
            fout.write(
                "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\t"
                "marker\tevalue\tbitscore\tuscg\ttaxid\ttaxname\tannotation_accessions\t"
                "annotation_description\n"
            )
            for line in utils.read_file(
                outputs.annotate_genes_output, skip_header=True
            ):
                contig = line.strip().split("\t")[0].rsplit("_", 1)[0]
                if contig in provirus_dict:
                    start = int(line.strip().split("\t")[1])
                    end = int(line.strip().split("\t")[2])
                    in_provirus = np.array(
                        [
                            start >= p.start and end <= p.end
                            for p in provirus_dict[contig]
                        ]
                    )
                    if any(in_provirus):
                        provirus_name = provirus_dict[contig][
                            np.argwhere(in_provirus)[0][0]
                        ].provirus_name
                        gene_number = line.strip().split("\t")[0].rsplit("_", 1)[1]
                        gene_name = f"{provirus_name}_{gene_number}"
                        fout.write(f"{gene_name}\t")
                        fout.write("\t".join(line.strip().split("\t")[1:]))
                        fout.write("\n")
        console.log(
            f"Provirus gene data was written to [green]{outputs.find_proviruses_genes_output.name}[/green]."
        )

    # Write `find_proviruses_taxonomy_output`
    with console.status(
        f"Writing [green]{outputs.find_proviruses_taxonomy_output.name}[/green]."
    ):
        taxonomy.write_taxonomic_assignment(
            outputs.find_proviruses_taxonomy_output,
            outputs.find_proviruses_genes_output,
            database_obj,
        )
        console.log(
            f"Taxonomic assignment data was written to [green]{outputs.find_proviruses_taxonomy_output.name}[/green]."
        )

    console.log("geNomad find-proviruses finished!", style="yellow")
