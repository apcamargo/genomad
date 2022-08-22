import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import numpy as np
import xgboost as xgb
from genomad import database, sequence, utils
from genomad._paths import GenomadData, GenomadOutputs

# Ignore numpy's overflow warning. For very long sequences, the `v_vs_c_score`,
# `v_vs_p_score` and `p_vs_c_score` values can be very high before the logistic
# transformation, causing an overflow in the np.exp function
np.warnings.filterwarnings("ignore", "overflow")

FEATURE_FILE_HEADER = "\t".join(
    [
        "seq_name",
        "n_genes",
        "n_uscg",
        "n_plasmid_hallmarks",
        "n_virus_hallmarks",
        "genetic_code",
        "strand_switch_rate",
        "coding_density",
        "no_rbs_freq",
        "sd_bacteroidetes_rbs_freq",
        "sd_canonical_rbs_freq",
        "tatata_rbs_freq",
        "cc_marker_freq",
        "cp_marker_freq",
        "cv_marker_freq",
        "pc_marker_freq",
        "pp_marker_freq",
        "pv_marker_freq",
        "vc_marker_freq",
        "vp_marker_freq",
        "vv_marker_freq",
        "c_marker_freq",
        "p_marker_freq",
        "v_marker_freq",
        "median_c_spm",
        "median_p_spm",
        "median_v_spm",
        "v_vs_c_score_logistic",
        "v_vs_p_score_logistic",
        "p_vs_c_score_logistic",
        "gv_marker_freq",
        "marker_enrichment_c",
        "marker_enrichment_p",
        "marker_enrichment_v",
    ]
)


@dataclass
class AnnotatedContig:
    seq_name: str
    contig_length: int
    coding_length: int = 0
    n_genes: int = 0
    n_usgc: int = 0
    n_plasmid_hallmarks: int = 0
    n_virus_hallmarks: int = 0
    genetic_code: int = 11
    n_cc_markers: int = 0
    n_cp_markers: int = 0
    n_cv_markers: int = 0
    n_pc_markers: int = 0
    n_pp_markers: int = 0
    n_pv_markers: int = 0
    n_vc_markers: int = 0
    n_vp_markers: int = 0
    n_vv_markers: int = 0
    n_gv_markers: int = 0
    spm_c: List[float] = field(default_factory=list)
    spm_p: List[float] = field(default_factory=list)
    spm_v: List[float] = field(default_factory=list)
    gene_strands: List[int] = field(default_factory=list)
    gene_rbs: List[str] = field(default_factory=list)

    @property
    def n_c_markers(self) -> int:
        return self.n_cc_markers + self.n_cp_markers + self.n_cv_markers

    @property
    def n_p_markers(self) -> int:
        return self.n_pc_markers + self.n_pp_markers + self.n_pv_markers

    @property
    def n_v_markers(self) -> int:
        return self.n_vc_markers + self.n_vp_markers + self.n_vv_markers

    @property
    def n_markers(self) -> int:
        return self.n_c_markers + self.n_p_markers + self.n_v_markers

    @property
    def coding_density(self) -> int:
        return self.coding_length / self.contig_length

    @property
    def strand_switch_rate(self) -> float:
        switches = sum(
            self.gene_strands[i] != self.gene_strands[i + 1]
            for i in range(self.n_genes - 1)
        )
        return switches / (self.n_genes - 1) if self.n_genes >= 2 else 0.0

    def get_compound_score(self, score_type: str) -> float:
        possible_score_types = {"v_vs_c", "v_vs_p", "p_vs_c"}
        if score_type not in possible_score_types:
            raise ValueError(
                f"Invalid score type. Expected one of: {possible_score_types}"
            )
        if score_type == "v_vs_c":
            score = (np.exp(self.spm_v) - np.exp(self.spm_c)).sum()
        elif score_type == "v_vs_p":
            score = (np.exp(self.spm_v) - np.exp(self.spm_p)).sum()
        elif score_type == "p_vs_c":
            score = (np.exp(self.spm_p) - np.exp(self.spm_c)).sum()
        return score

    def get_marker_enrichment(self, score_type: str) -> float:
        possible_score_types = {"c", "p", "v"}
        if score_type not in possible_score_types:
            raise ValueError(
                f"Invalid score type. Expected one of: {possible_score_types}"
            )
        if score_type == "c":
            score = (np.exp(self.spm_c) - np.exp(np.add(self.spm_p, self.spm_v))).sum()
        elif score_type == "p":
            score = (np.exp(self.spm_p) - np.exp(np.add(self.spm_c, self.spm_v))).sum()
        elif score_type == "v":
            score = (np.exp(self.spm_v) - np.exp(np.add(self.spm_c, self.spm_p))).sum()
        return score


def yield_annotated_contigs(
    input_path: Path,
    genes_output: Path,
    database_obj: database.Database,
    rbs_categories_dict: dict,
) -> AnnotatedContig:
    annotated_contigs_dict = {
        seq.accession: AnnotatedContig(seq.accession, len(seq))
        for seq in sequence.read_fasta(input_path, strip_n=True)
    }
    marker_features_dict = database_obj.get_marker_features()
    for line in utils.read_file(genes_output, skip_header=True):
        (
            gene,
            _,
            _,
            gene_length,
            strand,
            _,
            genetic_code,
            rbs,
            match,
            *_,
        ) = line.strip().split("\t")
        contig, gene_length, strand, genetic_code = (
            gene.rsplit("_", 1)[0],
            int(gene_length),
            int(strand),
            int(genetic_code),
        )
        (
            specificity_class,
            spm_c,
            spm_p,
            spm_v,
            gv_marker,
            uscg,
            plasmid_hallmark,
            virus_hallmark,
        ) = marker_features_dict.get(match, (None, 0.0, 0.0, 0.0, 0, 0, 0, 0))
        # Contigs that only contain Ns won't be in `annotated_contigs_dict`
        if contig in annotated_contigs_dict:
            annotated_contigs_dict[contig].n_genes += 1
            annotated_contigs_dict[contig].coding_length += gene_length
            annotated_contigs_dict[contig].gene_strands.append(strand)
            annotated_contigs_dict[contig].gene_rbs.append(rbs_categories_dict[rbs])
            annotated_contigs_dict[contig].genetic_code = genetic_code
            if specificity_class:
                annotated_contigs_dict[contig].spm_c.append(spm_c)
                annotated_contigs_dict[contig].spm_p.append(spm_p)
                annotated_contigs_dict[contig].spm_v.append(spm_v)
                annotated_contigs_dict[contig].n_gv_markers += gv_marker
                annotated_contigs_dict[contig].n_usgc += uscg
                annotated_contigs_dict[contig].n_plasmid_hallmarks += plasmid_hallmark
                annotated_contigs_dict[contig].n_virus_hallmarks += virus_hallmark
                if specificity_class == "CC":
                    annotated_contigs_dict[contig].n_cc_markers += 1
                elif specificity_class == "CP":
                    annotated_contigs_dict[contig].n_cp_markers += 1
                elif specificity_class == "CV":
                    annotated_contigs_dict[contig].n_cv_markers += 1
                elif specificity_class == "PC":
                    annotated_contigs_dict[contig].n_pc_markers += 1
                elif specificity_class == "PP":
                    annotated_contigs_dict[contig].n_pp_markers += 1
                elif specificity_class == "PV":
                    annotated_contigs_dict[contig].n_pv_markers += 1
                elif specificity_class == "VC":
                    annotated_contigs_dict[contig].n_vc_markers += 1
                elif specificity_class == "VP":
                    annotated_contigs_dict[contig].n_vp_markers += 1
                elif specificity_class == "VV":
                    annotated_contigs_dict[contig].n_vv_markers += 1
    yield from annotated_contigs_dict.values()


def get_feature_array(
    input_path: Path,
    genes_output: Path,
    database_obj: database.Database,
    rbs_file: Path,
):
    """
    The features in the array are, in order:
    (1) strand switch rate; (2) coding density; (3) no RBS freq.;
    (4) SD Bacteroidetes RBS freq.; (5) SD Canonical RBS freq.;
    (6) TATATA RBS frq.; (7) CC freq.; (8) CP freq.; (9) CV freq.;
    (10) PC freq.; (11) PP freq.; (12) PV freq.; (13) VC freq.;
    (14) VP freq.; (15) VV freq.; (16) C freq.; (17) P freq.;
    (18) V freq.; (19) median C SPM; (20) median P SPM;
    (21) median V SPM; (22) V vs. C score; (23) V vs. P score;
    (24) P vs. C score; (25) GV marker freq.
    """
    contig_array = []
    n_genes_array = []
    n_uscg_array = []
    n_hallmarks_array = []
    genetic_code_array = []
    feature_array = []
    marker_enrichment_array = []

    rbs_categories_dict = {}
    for line in utils.read_file(rbs_file):
        rbs, category = line.strip().split("\t")
        rbs_categories_dict[rbs] = category

    for annotated_contig in yield_annotated_contigs(
        input_path, genes_output, database_obj, rbs_categories_dict
    ):
        contig_array.append(annotated_contig.seq_name)
        rbs_frequency_dict = {
            k: v / annotated_contig.n_genes
            for k, v in Counter(annotated_contig.gene_rbs).items()
        }
        n_genes_array.append(annotated_contig.n_genes)
        n_uscg_array.append(annotated_contig.n_usgc)
        n_hallmarks_array.append([annotated_contig.n_plasmid_hallmarks, annotated_contig.n_virus_hallmarks])
        genetic_code_array.append(annotated_contig.genetic_code)
        contig_features = [
            annotated_contig.strand_switch_rate,
            annotated_contig.coding_density,
            rbs_frequency_dict.get("None", 0.0),
            rbs_frequency_dict.get("SD_Bacteroidetes", 0.0),
            rbs_frequency_dict.get("SD_Canonical", 0.0),
            rbs_frequency_dict.get("TATATA_3_6", 0.0),
            annotated_contig.n_cc_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_cp_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_cv_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_pc_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_pp_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_pv_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_vc_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_vp_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_vv_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_c_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_p_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            annotated_contig.n_v_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
            np.median(annotated_contig.spm_c)
            if annotated_contig.n_markers > 0
            else 0.0,
            np.median(annotated_contig.spm_p)
            if annotated_contig.n_markers > 0
            else 0.0,
            np.median(annotated_contig.spm_v)
            if annotated_contig.n_markers > 0
            else 0.0,
            utils.logistic(annotated_contig.get_compound_score("v_vs_c"), 2),
            utils.logistic(annotated_contig.get_compound_score("v_vs_p"), 2),
            utils.logistic(annotated_contig.get_compound_score("p_vs_c"), 2),
            annotated_contig.n_gv_markers / annotated_contig.n_genes
            if annotated_contig.n_genes > 0
            else 0.0,
        ]
        feature_array.append(contig_features)
        marker_enrichment = [
            annotated_contig.get_marker_enrichment("c"),
            annotated_contig.get_marker_enrichment("p"),
            annotated_contig.get_marker_enrichment("v"),
        ]
        marker_enrichment_array.append(marker_enrichment)
    return (
        np.array(contig_array),
        np.array(n_genes_array),
        np.array(n_uscg_array),
        np.array(n_hallmarks_array),
        np.array(genetic_code_array),
        np.array(feature_array),
        np.array(marker_enrichment_array),
    )


def main(input_path, output_path, database_path, restart, threads, verbose):
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
        output_file=outputs.marker_classification_log, verbose=verbose
    )
    # Create a dictionary containing the parameters
    parameter_dict = {}

    # Check if the find-provirus module was executed
    classify_proviruses = utils.check_provirus_execution(
        prefix, input_path, output_path
    )

    # Display the module header
    output_files_path_list = [
        outputs.marker_classification_execution_info,
        outputs.features_output,
        outputs.features_npz_output,
        outputs.marker_classification_output,
        outputs.marker_classification_npz_output,
    ]
    output_files_description_list = [
        "execution parameters",
        "sequence feature data: tabular format",
        "sequence feature data: binary format",
        "sequence classification: tabular format",
        "sequence classification: binary format",
    ]
    if classify_proviruses:
        output_files_path_list.extend(
            [
                outputs.provirus_features_output,
                outputs.provirus_features_npz_output,
                outputs.provirus_marker_classification_output,
                outputs.provirus_marker_classification_npz_output,
            ]
        )
        output_files_description_list.extend(
            [
                "provirus feature data: tabular format",
                "provirus feature data: binary format",
                "provirus classification: tabular format",
                "provirus classification: binary format",
            ]
        )

    utils.display_header(
        console,
        "marker-classification",
        (
            "This will classify the input sequences into chromosome, plasmid, "
            "or virus based on the presence of geNomad markers and other "
            "gene-related features."
        ),
        outputs.marker_classification_dir,
        output_files_path_list,
        output_files_description_list,
    )

    # Check if the annotation module was previously executed
    if not outputs.annotate_genes_output.exists():
        console.error(
            f"[u]{outputs.annotate_genes_output.name}[/u] was not found in the output "
            "directory. Please execute the [cyan]annotate[/cyan] module to generate it."
        )
        sys.exit(1)
    # Check if the input FASTA for `annotate` and `marker-classification` is the same
    elif not utils.compare_executions(
        input_path, {}, outputs.annotate_execution_info, only_md5=True
    ):
        console.error(
            "The input FASTA file is different from the one used in the [cyan]annotate[/cyan] "
            "module. Please execute both modules using the same input"
        )
        sys.exit(1)

    # Check the input FASTA file
    if not sequence.check_fasta(input_path):
        console.error(
            f"[u]{input_path}[/u] is either empty or contains multiple entries "
            "with the same identifier. Please check your input FASTA file and "
            "execute [cyan]genomad marker-classification[/cyan] again."
        )
        sys.exit(1)

    # Print initial log
    console.log("Executing [cyan]genomad marker-classification[/cyan].")

    # Check if steps can be skipped
    skip = False
    if (
        outputs.marker_classification_execution_info.exists()
        and any(path.exists() for path in output_files_path_list)
        and not restart
    ):
        # Check if the previous execution used the same input file and parameters
        if utils.compare_executions(
            input_path, parameter_dict, outputs.marker_classification_execution_info
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
    if not outputs.marker_classification_dir.is_dir():
        console.log(
            f"Creating the [green]{outputs.marker_classification_dir}[/green] directory."
        )
        outputs.marker_classification_dir.mkdir()

    # Write the execution data to `execution_info_output`
    utils.write_execution_info(
        "marker_classification",
        input_path,
        parameter_dict,
        outputs.marker_classification_execution_info,
    )

    # Initialize the database object
    database_obj = database.Database(database_path)

    # Compute features and write `features_npz_output`
    if skip and outputs.features_npz_output.exists():
        console.log(
            f"[green]{outputs.features_npz_output.name}[/green] was found. "
            "Skipping feature computation."
        )
        contig_names = np.load(outputs.features_npz_output)["contig_names"]
        contig_n_genes = np.load(outputs.features_npz_output)["contig_n_genes"]
        contig_n_uscg = np.load(outputs.features_npz_output)["contig_n_uscg"]
        contig_n_hallmarks = np.load(outputs.features_npz_output)["contig_n_hallmarks"]
        contig_genetic_code = np.load(outputs.features_npz_output)[
            "contig_genetic_code"
        ]
        contig_features = np.load(outputs.features_npz_output)["contig_features"]
        contig_marker_enrichment = np.load(outputs.features_npz_output)[
            "contig_marker_enrichment"
        ]
    else:
        with console.status("Computing sequence features."):
            (
                contig_names,
                contig_n_genes,
                contig_n_uscg,
                contig_n_hallmarks,
                contig_genetic_code,
                contig_features,
                contig_marker_enrichment,
            ) = get_feature_array(
                input_path,
                outputs.annotate_genes_output,
                database_obj,
                GenomadData.rbs_file,
            )
            console.log("Sequence features computed.")
        with console.status(
            "Writing sequence features in binary format to [green]"
            f"{outputs.features_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.features_npz_output,
                contig_names=contig_names,
                contig_n_genes=contig_n_genes,
                contig_n_uscg=contig_n_uscg,
                contig_n_hallmarks=contig_n_hallmarks,
                contig_genetic_code=contig_genetic_code,
                contig_features=contig_features,
                contig_marker_enrichment=contig_marker_enrichment,
            )
            console.log(
                "Sequence features in binary format written to [green]"
                f"{outputs.features_npz_output.name}[/green]."
            )

    # Write `features_output`
    with console.status(
        "Writing sequence features in tabular format to [green]"
        f"{outputs.features_output.name}[/green]."
    ):
        with open(outputs.features_output, "w") as fout:
            fout.write(f"{FEATURE_FILE_HEADER}\n")
            for (
                name,
                n_genes,
                n_uscg,
                n_hallmarks,
                genetic_code,
                features,
                marker_enrichment,
            ) in zip(
                contig_names,
                contig_n_genes,
                contig_n_uscg,
                contig_n_hallmarks,
                contig_genetic_code,
                contig_features,
                contig_marker_enrichment,
            ):
                features = "".join(map(lambda x: f"{x:.4f}\t", features)).strip()
                marker_enrichment = "".join(
                    map(lambda x: f"{x:.4f}\t", marker_enrichment)
                ).strip()
                fout.write(
                    f"{name}\t{n_genes}\t{n_uscg}\t{n_hallmarks[0]}\t{n_hallmarks[1]}\t"
                    f"{genetic_code}\t{features}\t{marker_enrichment}\n"
                )
        console.log(
            "Sequence features in tabular format written to [green]"
            f"{outputs.features_output.name}[/green]."
        )

    # Compute provirus features and write `provirus_features_npz_output`
    if skip and classify_proviruses and outputs.provirus_features_npz_output.exists():
        console.log(
            f"[green]{outputs.provirus_features_npz_output.name}[/green] was found. "
            "Skipping provirus feature computation."
        )
        provirus_names = np.load(outputs.provirus_features_npz_output)["provirus_names"]
        provirus_n_genes = np.load(outputs.provirus_features_npz_output)[
            "provirus_n_genes"
        ]
        provirus_n_uscg = np.load(outputs.provirus_features_npz_output)[
            "provirus_n_uscg"
        ]
        provirus_n_hallmarks = np.load(outputs.provirus_features_npz_output)["provirus_n_hallmarks"]
        provirus_genetic_code = np.load(outputs.provirus_features_npz_output)[
            "provirus_genetic_code"
        ]
        provirus_features = np.load(outputs.provirus_features_npz_output)[
            "provirus_features"
        ]
        provirus_marker_enrichment = np.load(outputs.provirus_features_npz_output)[
            "provirus_marker_enrichment"
        ]
    elif classify_proviruses:
        with console.status("Computing provirus features."):
            (
                provirus_names,
                provirus_n_genes,
                provirus_n_uscg,
                provirus_n_hallmarks,
                provirus_genetic_code,
                provirus_features,
                provirus_marker_enrichment,
            ) = get_feature_array(
                outputs.find_proviruses_nucleotide_output,
                outputs.find_proviruses_genes_output,
                database_obj,
                GenomadData.rbs_file,
            )
            console.log("Provirus features computed.")
        with console.status(
            "Writing provirus features in binary format to [green]"
            f"{outputs.provirus_features_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.provirus_features_npz_output,
                provirus_names=provirus_names,
                provirus_n_genes=provirus_n_genes,
                provirus_n_uscg=provirus_n_uscg,
                provirus_n_hallmarks=provirus_n_hallmarks,
                provirus_genetic_code=provirus_genetic_code,
                provirus_features=provirus_features,
                provirus_marker_enrichment=provirus_marker_enrichment,
            )
            console.log(
                "Provirus features in binary format written to [green]"
                f"{outputs.provirus_features_npz_output.name}[/green]."
            )

    # Write `provirus_features_output`
    if classify_proviruses:
        with console.status(
            "Writing provirus features in tabular format to [green]"
            f"{outputs.provirus_features_output.name}[/green]."
        ):
            with open(outputs.provirus_features_output, "w") as fout:
                fout.write(f"{FEATURE_FILE_HEADER}\n")
                for (
                    name,
                    n_genes,
                    n_uscg,
                    n_hallmarks,
                    genetic_code,
                    features,
                    marker_enrichment,
                ) in zip(
                    provirus_names,
                    provirus_n_genes,
                    provirus_n_uscg,
                    provirus_n_hallmarks,
                    provirus_genetic_code,
                    provirus_features,
                    provirus_marker_enrichment,
                ):
                    features = "".join(map(lambda x: f"{x:.4f}\t", features)).strip()
                    marker_enrichment = "".join(
                        map(lambda x: f"{x:.4f}\t", marker_enrichment)
                    ).strip()
                    fout.write(
                        f"{name}\t{n_genes}\t{n_uscg}\t{n_hallmarks[0]}\t{n_hallmarks[1]}\t"
                        f"{genetic_code}\t{features}\t{marker_enrichment}\n"
                    )
            console.log(
                "Provirus features in tabular format written to [green]"
                f"{outputs.provirus_features_output.name}[/green]."
            )

    # Classify sequences and write `marker_classification_npz_output`
    if skip and outputs.marker_classification_npz_output.exists():
        console.log(
            f"[green]{outputs.marker_classification_npz_output.name}[/green] was found. "
            "Skipping sequence classification."
        )
        contig_predictions = np.load(outputs.marker_classification_npz_output)[
            "predictions"
        ]
    else:
        if not len(contig_features):
            console.error("No sequences were found. Please check your input FASTA.")
            sys.exit(1)
        with console.status("Classifying sequences."):
            df_model = xgb.Booster(
                {"nthread": threads}, model_file=GenomadData.decision_forest_file
            )
            contig_predictions = utils.softmax(
                df_model.predict(xgb.DMatrix(contig_features), output_margin=True),
                temperature=2,
            )
            console.log("Sequences classified.")
        with console.status(
            "Writing sequence classification in binary format to [green]"
            f"{outputs.marker_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.marker_classification_npz_output,
                contig_names=contig_names,
                predictions=contig_predictions,
            )
            console.log(
                "Sequence classification in binary format written to "
                f"[green]{outputs.marker_classification_npz_output.name}[/green]."
            )

    # Write `marker_classification_output`
    with console.status(
        "Writing sequence classification in tabular format to [green]"
        f"{outputs.marker_classification_output.name}[/green]."
    ):
        with open(outputs.marker_classification_output, "w") as fout:
            fout.write("seq_name\tchromosome_score\tplasmid_score\tvirus_score\n")
            for name, scores in zip(contig_names, contig_predictions):
                scores = "".join(map(lambda x: f"{x:.4f}\t", scores)).strip()
                fout.write(f"{name}\t{scores}\n")
        console.log(
            "Sequence classification in tabular format written to [green]"
            f"{outputs.marker_classification_output.name}[/green]."
        )

    # Classify proviruses and write `provirus_marker_classification_npz_output`
    if (
        skip
        and classify_proviruses
        and outputs.provirus_marker_classification_npz_output.exists()
    ):
        console.log(
            f"[green]{outputs.provirus_marker_classification_npz_output.name}[/green] was found. "
            "Skipping provirus classification."
        )
        provirus_predictions = np.load(
            outputs.provirus_marker_classification_npz_output
        )["predictions"]
    elif classify_proviruses:
        with console.status("Classifying proviruses."):
            df_model = xgb.Booster(
                {"nthread": threads}, model_file=GenomadData.decision_forest_file
            )
            provirus_predictions = utils.softmax(
                df_model.predict(xgb.DMatrix(provirus_features), output_margin=True),
                temperature=2,
            )
            console.log("Proviruses classified.")
        with console.status(
            "Writing provirus classification in binary format to [green]"
            f"{outputs.provirus_marker_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.provirus_marker_classification_npz_output,
                provirus_names=provirus_names,
                predictions=provirus_predictions,
            )
            console.log(
                "Provirus classification in binary format written to "
                f"[green]{outputs.provirus_marker_classification_npz_output.name}[/green]."
            )

    # Write `provirus_marker_classification_output`
    if classify_proviruses:
        with console.status(
            "Writing provirus classification in tabular format to [green]"
            f"{outputs.provirus_marker_classification_output.name}[/green]."
        ):
            with open(outputs.provirus_marker_classification_output, "w") as fout:
                fout.write("seq_name\tchromosome_score\tplasmid_score\tvirus_score\n")
                for name, scores in zip(provirus_names, provirus_predictions):
                    scores = "".join(map(lambda x: f"{x:.4f}\t", scores)).strip()
                    fout.write(f"{name}\t{scores}\n")
            console.log(
                "Provirus classification in tabular format written to [green]"
                f"{outputs.provirus_marker_classification_output.name}[/green]."
            )

    console.log("geNomad marker-classification finished!", style="yellow")
