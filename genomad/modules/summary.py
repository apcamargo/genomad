import itertools
import sys
from collections import defaultdict

import numpy as np
from genomad import sequence, utils
from genomad._paths import GenomadOutputs


def get_fdr_array(probability_array):
    n_errors = 0
    fdr_array = []
    for n, p in enumerate(probability_array, start=1):
        p_error = 1 - p
        n_errors += p_error
        fdr_array.append(n_errors / n)
    return np.array(fdr_array)


def flag_sequences(
    contig_name_array,
    contig_score_array,
    class_index,
    min_score,
    max_fdr,
    min_marker_enrichment,
    min_hallmarks,
    max_uscg,
    filters_dict,
    provirus_name_array=None,
    provirus_score_array=None,
):
    if (
        (provirus_name_array is not None)
        and (provirus_score_array is not None)
        and (len(provirus_name_array) and len(provirus_score_array))
    ):
        name_array = np.concatenate([contig_name_array, provirus_name_array])
        score_array = np.concatenate([contig_score_array, provirus_score_array])
        provirus_name_set = set(provirus_name_array)
    else:
        name_array = contig_name_array
        score_array = contig_score_array
        provirus_name_set = set()
    selected_name_array = []
    selected_score_array = []
    added_contigs = set()
    added_proviruses = set()
    for i in reversed(score_array[:, class_index].argsort()):
        n_uscg, marker_enrichment, n_hallmarks = filters_dict.get(
            name_array[i], (0, np.zeros(3), (0, 0))
        )
        marker_enrichment = marker_enrichment[class_index]
        n_hallmarks = n_hallmarks[class_index - 1]
        if (
            score_array[i].argmax() == class_index
            and score_array[i, class_index] >= min_score
            and marker_enrichment >= min_marker_enrichment
            and n_hallmarks >= min_hallmarks
            and n_uscg <= max_uscg
        ):
            if name_array[i] in provirus_name_set:
                contig_name = name_array[i].rsplit("|", 1)[0]
                if contig_name not in added_contigs:
                    selected_name_array.append(name_array[i])
                    selected_score_array.append(score_array[i, class_index])
                    added_proviruses.add(contig_name)
            else:
                contig_name = name_array[i]
                if contig_name not in added_proviruses:
                    selected_name_array.append(name_array[i])
                    selected_score_array.append(score_array[i, class_index])
                    added_contigs.add(contig_name)
    if not max_fdr:
        return (
            np.array(selected_name_array),
            np.array(selected_score_array),
            np.array([]),
        )
    fdr_array = get_fdr_array(selected_score_array)
    return (
        np.array(selected_name_array)[fdr_array < max_fdr],
        np.array(selected_score_array)[fdr_array < max_fdr],
        fdr_array[fdr_array < max_fdr],
    )


def main(
    input_path,
    output_path,
    verbose,
    min_score,
    max_fdr,
    min_plasmid_marker_enrichment,
    min_virus_marker_enrichment,
    min_plasmid_hallmarks,
    min_virus_hallmarks,
    max_uscg,
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
    console = utils.HybridConsole(output_file=outputs.summary_log, verbose=verbose)
    # Create a dictionary containing the parameters
    parameter_dict = {
        "min_score": min_score,
        "max_fdr": max_fdr,
        "min_plasmid_hallmarks": min_plasmid_hallmarks,
        "min_virus_hallmarks": min_virus_hallmarks,
        "min_plasmid_marker_enrichment": min_plasmid_marker_enrichment,
        "min_virus_marker_enrichment": min_virus_marker_enrichment,
        "max_uscg": max_uscg,
    }

    # Check which modules were previously executed
    annotate_exec = bool(
        outputs.annotate_execution_info.exists()
        and outputs.annotate_proteins_output.exists()
        and outputs.annotate_genes_output.exists()
        and outputs.annotate_taxonomy_output.exists()
    )
    marker_classification_exec = bool(
        outputs.marker_classification_execution_info.exists()
        and outputs.marker_classification_npz_output.exists()
        and outputs.features_npz_output.exists()
    )
    nn_classification_exec = bool(
        outputs.nn_classification_execution_info.exists()
        and outputs.nn_classification_npz_output.exists()
    )
    aggregated_classification_exec = bool(
        outputs.aggregated_classification_execution_info.exists()
        and outputs.aggregated_classification_npz_output.exists()
    )
    find_proviruses_exec = bool(
        outputs.find_proviruses_execution_info.exists()
        and outputs.find_proviruses_output.exists()
        and outputs.find_proviruses_nucleotide_output.exists()
        and outputs.find_proviruses_proteins_output.exists()
        and outputs.find_proviruses_genes_output.exists()
    )
    marker_classification_provirus_exec = bool(
        marker_classification_exec
        and find_proviruses_exec
        and outputs.provirus_marker_classification_npz_output.exists()
    )
    nn_classification_provirus_exec = bool(
        nn_classification_exec
        and find_proviruses_exec
        and outputs.provirus_nn_classification_npz_output.exists()
    )
    aggregated_classification_provirus_exec = bool(
        aggregated_classification_exec
        and find_proviruses_exec
        and outputs.provirus_aggregated_classification_npz_output.exists()
    )
    score_calibration_marker_exec = bool(
        outputs.score_calibration_execution_info.exists()
        and outputs.calibrated_marker_classification_output.exists()
    )
    score_calibration_nn_exec = bool(
        outputs.score_calibration_execution_info.exists()
        and outputs.calibrated_nn_classification_npz_output.exists()
    )
    score_calibration_aggregated_exec = bool(
        outputs.score_calibration_execution_info.exists()
        and outputs.calibrated_aggregated_classification_npz_output.exists()
    )
    score_calibration_marker_provirus_exec = bool(
        outputs.score_calibration_execution_info.exists()
        and outputs.provirus_calibrated_marker_classification_npz_output.exists()
    )
    score_calibration_nn_provirus_exec = bool(
        outputs.score_calibration_execution_info.exists()
        and outputs.provirus_calibrated_nn_classification_npz_output.exists()
    )
    score_calibration_aggregated_provirus_exec = bool(
        outputs.score_calibration_execution_info.exists()
        and outputs.provirus_calibrated_aggregated_classification_npz_output.exists()
    )

    # Pick the classification that will be used in the summary files
    if score_calibration_aggregated_exec:
        selected_classifier = "calibrated_aggregated"
        include_provirus = True if score_calibration_aggregated_provirus_exec else False
        contig_classification_output = (
            outputs.calibrated_aggregated_classification_npz_output
        )
        provirus_classification_output = (
            outputs.provirus_calibrated_aggregated_classification_npz_output
        )
    elif aggregated_classification_exec:
        selected_classifier = "aggregated"
        include_provirus = True if aggregated_classification_provirus_exec else False
        contig_classification_output = outputs.aggregated_classification_npz_output
        provirus_classification_output = (
            outputs.provirus_aggregated_classification_npz_output
        )
    elif score_calibration_marker_exec:
        selected_classifier = "calibrated_marker"
        include_provirus = True if score_calibration_marker_provirus_exec else False
        contig_classification_output = (
            outputs.calibrated_marker_classification_npz_output
        )
        provirus_classification_output = (
            outputs.provirus_calibrated_marker_classification_npz_output
        )
    elif marker_classification_exec:
        selected_classifier = "marker"
        include_provirus = True if marker_classification_provirus_exec else False
        contig_classification_output = outputs.marker_classification_npz_output
        provirus_classification_output = (
            outputs.provirus_marker_classification_npz_output
        )
    elif score_calibration_nn_exec:
        selected_classifier = "calibrated_nn"
        include_provirus = True if score_calibration_nn_provirus_exec else False
        contig_classification_output = outputs.calibrated_nn_classification_npz_output
        provirus_classification_output = (
            outputs.provirus_calibrated_nn_classification_npz_output
        )
    elif nn_classification_exec:
        selected_classifier = "nn"
        include_provirus = True if nn_classification_provirus_exec else False
        contig_classification_output = outputs.nn_classification_npz_output
        provirus_classification_output = (
            outputs.provirus_calibrated_nn_classification_npz_output
        )
    else:
        console.error(
            "No previous execution of the [cyan]marker-classification[/cyan], "
            "[cyan]nn-classification[/cyan], [cyan]aggregated-classification[/cyan], "
            "or [cyan]score-calibration[/cyan] were detected. Please execute at "
            "least one of these modules."
        )
        sys.exit(1)

    # Display header
    output_files_path_list = [
        outputs.summary_execution_info,
        outputs.summary_virus_output,
        outputs.summary_plasmid_output,
        outputs.summary_virus_sequences_output,
        outputs.summary_plasmid_sequences_output,
    ]
    output_files_description_list = [
        "execution parameters",
        "virus classification summary",
        "plasmid classification summary",
        "virus nucleotide FASTA file",
        "plasmid nucleotide FASTA file",
    ]
    if annotate_exec:
        output_files_path_list.extend(
            [
                outputs.summary_virus_proteins_output,
                outputs.summary_plasmid_proteins_output,
                outputs.summary_virus_genes_output,
                outputs.summary_plasmid_genes_output,
            ]
        )
        output_files_description_list.extend(
            [
                "virus protein FASTA file",
                "plasmid protein FASTA file",
                "virus gene annotation data",
                "plasmid gene annotation data",
            ]
        )
    utils.display_header(
        console,
        "summary",
        (
            "This will summarize the results across modules into a classification report."
        ),
        outputs.summary_dir,
        output_files_path_list,
        output_files_description_list,
    )

    # Check if the inputs of the previous modules were all the same
    md5_list = [utils.get_md5(input_path)]
    if find_proviruses_exec:
        md5_list.append(
            utils.get_execution_info(outputs.find_proviruses_execution_info)[0]
        )
    if marker_classification_exec:
        md5_list.append(
            utils.get_execution_info(outputs.marker_classification_execution_info)[0]
        )
    if nn_classification_exec:
        md5_list.append(
            utils.get_execution_info(outputs.nn_classification_execution_info)[0]
        )
    if aggregated_classification_exec:
        md5_list.append(
            utils.get_execution_info(outputs.aggregated_classification_execution_info)[
                0
            ]
        )
    if (
        score_calibration_marker_exec
        or score_calibration_nn_exec
        or score_calibration_aggregated_exec
        or score_calibration_marker_provirus_exec
        or score_calibration_nn_provirus_exec
        or score_calibration_aggregated_provirus_exec
    ):
        md5_list.append(
            utils.get_execution_info(outputs.score_calibration_execution_info)[0]
        )
    if len(set(md5_list)) > 1:
        console.error(
            "Different input FASTA files were used as input for the different modules. Please "
            "execute all the modules using the same input."
        )
        sys.exit(1)

    # Print initial log
    console.log("Executing [cyan]genomad summary[/cyan].")

    # If necessary, create the subdirectory
    if not outputs.summary_dir.is_dir():
        console.log(f"Creating the [green]{outputs.summary_dir}[/green] directory.")
        outputs.summary_dir.mkdir()

    # Write the execution data to `execution_info_output`
    utils.write_execution_info(
        "summary",
        input_path,
        parameter_dict,
        outputs.summary_execution_info,
    )

    if selected_classifier == "calibrated_aggregated":
        console.log(
            "Using calibrated scores from [cyan]aggregated-classification[/cyan]."
        )
    elif selected_classifier == "aggregated":
        console.log("Using scores from [cyan]aggregated-classification[/cyan].")
    elif selected_classifier == "calibrated_marker":
        console.log("Using calibrated scores from [cyan]marker-classification[/cyan].")
    elif selected_classifier == "aggregated":
        console.log("Using scores from [cyan]marker-classification[/cyan].")
    elif selected_classifier == "calibrated_nn":
        console.log("Using calibrated scores from [cyan]nn-classification[/cyan].")
    elif selected_classifier == "nn":
        console.log("Using scores from [cyan]nn-classification[/cyan].")

    # Warn if aggregated-classification was not executed
    if (
        (marker_classification_exec and nn_classification_exec)
        and not aggregated_classification_exec
    ) or (
        (marker_classification_exec and score_calibration_marker_exec)
        and not score_calibration_aggregated_exec
    ):
        console.warning(
            "The outputs of [cyan]marker-classification[/cyan] and "
            "[cyan]nn-classification[/cyan] were detected, but "
            "[cyan]aggregated-classification[/cyan] was not executed. The classification "
            "summary will prioritize [cyan]marker-classification[/cyan], but "
            "it is recommended to execute [cyan]aggregated-classification[/cyan] before "
            "[cyan]summary[/cyan] to improve the classification performance."
        )

    # Store the number of USCGs, genetic code, and number of genes
    n_genes_dict = {}
    genetic_code_dict = {}
    filters_dict = {}
    if marker_classification_exec:
        for k, v1, v2, v3, v4, v5 in zip(
            np.load(outputs.features_npz_output)["contig_names"],
            np.load(outputs.features_npz_output)["contig_n_uscg"],
            np.load(outputs.features_npz_output)["contig_n_genes"],
            np.load(outputs.features_npz_output)["contig_genetic_code"],
            np.load(outputs.features_npz_output)["contig_marker_enrichment"],
            np.load(outputs.features_npz_output)["contig_n_hallmarks"],
        ):
            n_genes_dict[k] = v2
            genetic_code_dict[k] = v3
            filters_dict[k] = (v1, v4, v5)
        if include_provirus:
            for k, v1, v2, v3, v4, v5 in zip(
                np.load(outputs.provirus_features_npz_output)["provirus_names"],
                np.load(outputs.provirus_features_npz_output)["provirus_n_uscg"],
                np.load(outputs.provirus_features_npz_output)["provirus_n_genes"],
                np.load(outputs.provirus_features_npz_output)["provirus_genetic_code"],
                np.load(outputs.provirus_features_npz_output)[
                    "provirus_marker_enrichment"
                ],
                np.load(outputs.provirus_features_npz_output)["provirus_n_hallmarks"],
            ):
                n_genes_dict[k] = v2
                genetic_code_dict[k] = v3
                filters_dict[k] = (v1, v4, v5)

    # Load the classificaion results:
    score_dict = {}
    score_dict["contig_names"] = np.load(contig_classification_output)["contig_names"]
    score_dict["contig_predictions"] = np.load(contig_classification_output)[
        "predictions"
    ]
    if include_provirus:
        score_dict["provirus_names"] = np.load(provirus_classification_output)[
            "provirus_names"
        ]
        score_dict["provirus_predictions"] = np.load(provirus_classification_output)[
            "predictions"
        ]
    else:
        score_dict["provirus_names"] = np.array([])
        score_dict["provirus_predictions"] = np.array([])

    # Flag viruses and plasmids
    with console.status("Identifying plasmids and viruses."):
        if not selected_classifier.startswith("calibrated"):
            max_fdr = None
        plasmid_name_array, plasmid_score_array, plasmid_fdr_array = flag_sequences(
            score_dict["contig_names"],
            score_dict["contig_predictions"],
            1,
            min_score,
            max_fdr,
            min_plasmid_marker_enrichment,
            min_plasmid_hallmarks,
            max_uscg,
            filters_dict,
        )
        virus_name_array, virus_score_array, virus_fdr_array = flag_sequences(
            score_dict["contig_names"],
            score_dict["contig_predictions"],
            2,
            min_score,
            max_fdr,
            min_virus_marker_enrichment,
            min_virus_hallmarks,
            max_uscg,
            filters_dict,
            score_dict["provirus_names"],
            score_dict["provirus_predictions"],
        )
        plasmid_name_set = set(plasmid_name_array)
        virus_name_set = set(virus_name_array)
        console.log(
            f"{len(plasmid_name_array):,} plasmid(s) and "
            f"{len(virus_name_array):,} virus(es) were identified."
        )

    with console.status("Writing nucleotide FASTA files."):
        length_dict = {}
        terminal_repeat_dict = {}
        with open(outputs.summary_plasmid_sequences_output, "w") as fout1, open(
            outputs.summary_virus_sequences_output, "w"
        ) as fout2:
            for seq in sequence.read_fasta(input_path):
                if seq.accession in plasmid_name_set:
                    length_dict[seq.accession] = len(seq)
                    if seq.has_dtr():
                        terminal_repeat_dict[seq.accession] = "DTR"
                    elif seq.has_itr():
                        terminal_repeat_dict[seq.accession] = "ITR"
                    else:
                        terminal_repeat_dict[seq.accession] = "No terminal repeats"
                    fout1.write(str(seq))
                elif seq.accession in virus_name_set:
                    length_dict[seq.accession] = len(seq)
                    if seq.has_dtr():
                        terminal_repeat_dict[seq.accession] = "DTR"
                    elif seq.has_itr():
                        terminal_repeat_dict[seq.accession] = "ITR"
                    else:
                        terminal_repeat_dict[seq.accession] = "No terminal repeats"
                    fout2.write(str(seq))
            if include_provirus:
                for seq in sequence.read_fasta(
                    outputs.find_proviruses_nucleotide_output
                ):
                    if seq.accession in virus_name_set:
                        length_dict[seq.accession] = len(seq)
                        terminal_repeat_dict[seq.accession] = "Provirus"
                        fout2.write(str(seq))
        console.log(
            "Nucleotide sequences were written to "
            f"[green]{outputs.summary_plasmid_sequences_output.name}[/green] and "
            f"[green]{outputs.summary_virus_sequences_output.name}[/green]."
        )

    if annotate_exec:
        with console.status("Writing protein FASTA files."):
            with open(outputs.summary_plasmid_proteins_output, "w") as fout1, open(
                outputs.summary_virus_proteins_output, "w"
            ) as fout2:
                for seq in sequence.read_fasta(outputs.annotate_proteins_output):
                    if seq.accession.rsplit("_", 1)[0] in plasmid_name_set:
                        fout1.write(str(seq))
                    elif seq.accession.rsplit("_", 1)[0] in virus_name_set:
                        fout2.write(str(seq))
                if include_provirus:
                    for seq in sequence.read_fasta(
                        outputs.find_proviruses_proteins_output
                    ):
                        if seq.accession.rsplit("_", 1)[0] in virus_name_set:
                            fout2.write(str(seq))
            console.log(
                "Protein sequences were written to "
                f"[green]{outputs.summary_plasmid_proteins_output.name}[/green] and "
                f"[green]{outputs.summary_virus_proteins_output.name}[/green]."
            )

        with console.status("Writting gene annotation data."):
            conjscan_genes_dict = defaultdict(list)
            amr_genes_dict = defaultdict(list)
            with open(outputs.summary_plasmid_genes_output, "w") as fout1, open(
                outputs.summary_virus_genes_output, "w"
            ) as fout2:
                fout1.write(
                    "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\tmarker\t"
                    "evalue\tbitscore\tuscg\tplasmid_hallmark\tvirus_hallmark\ttaxid\ttaxname\t"
                    "annotation_conjscan\tannotation_amr\tannotation_accessions\tannotation_description\n"
                )
                fout2.write(
                    "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\tmarker\t"
                    "evalue\tbitscore\tuscg\tplasmid_hallmark\tvirus_hallmark\ttaxid\ttaxname\t"
                    "annotation_conjscan\tannotation_amr\tannotation_accessions\tannotation_description\n"
                )
                for line in utils.read_file(
                    outputs.annotate_genes_output, skip_header=True
                ):
                    seq_name = line.split("\t")[0].rsplit("_", 1)[0]
                    if seq_name in plasmid_name_set:
                        fout1.write(line)
                        if line.split("\t")[16] != "NA":
                            conjscan_genes_dict[seq_name].append(line.split("\t")[16])
                        if line.split("\t")[17] != "NA":
                            amr_genes_dict[seq_name].append(line.split("\t")[17])
                    elif seq_name in virus_name_set:
                        fout2.write(line)
                if include_provirus:
                    for line in utils.read_file(
                        outputs.find_proviruses_genes_output, skip_header=True
                    ):
                        if line.split("\t")[0].rsplit("_", 1)[0] in virus_name_set:
                            fout2.write(line)
            console.log(
                "Gene annotation data was written to "
                f"[green]{outputs.summary_plasmid_genes_output.name}[/green] and "
                f"[green]{outputs.summary_virus_genes_output.name}[/green]."
            )

    with console.status("Writing summary files."):
        # Store provirus coordinates
        provirus_coord_dict = {}
        if include_provirus:
            for line in utils.read_file(
                outputs.find_proviruses_output, skip_header=True
            ):
                seq_name, _, start, end, *_ = line.strip().split("\t")
                if seq_name in virus_name_set:
                    start, end = int(start), int(end)
                    provirus_coord_dict[seq_name] = (start, end)
        # Store the predicted taxonomy
        taxonomy_dict = {}
        if annotate_exec:
            for line in utils.read_file(
                outputs.annotate_taxonomy_output, skip_header=True
            ):
                seq_name, _, _, _, lineage = line.strip().split("\t")
                if seq_name in virus_name_set:
                    taxonomy_dict[seq_name] = lineage
            if include_provirus:
                for line in utils.read_file(
                    outputs.find_proviruses_taxonomy_output, skip_header=True
                ):
                    seq_name, _, _, _, lineage = line.strip().split("\t")
                    if seq_name in virus_name_set:
                        taxonomy_dict[seq_name] = lineage

        # Write plasmid summary file
        with open(outputs.summary_plasmid_output, "w") as fout:
            fout.write(
                "seq_name\tlength\ttopology\tn_genes\tgenetic_code\tplasmid_score\t"
                "fdr\tn_hallmarks\tmarker_enrichment\tconjugation_genes\tamr_genes\n"
            )
            for seq_name, score, fdr in itertools.zip_longest(
                plasmid_name_array,
                plasmid_score_array,
                plasmid_fdr_array,
                fillvalue="NA",
            ):
                length = length_dict.get(seq_name, "NA")
                terminal_repeat = terminal_repeat_dict.get(seq_name, "NA")
                n_genes = n_genes_dict.get(seq_name, "NA")
                genetic_code = genetic_code_dict.get(seq_name, "NA")
                score = f"{score:.4f}"
                fdr = f"{fdr:.4f}" if not isinstance(fdr, str) else fdr
                if annotate_exec:
                    _, marker_enrichment, n_hallmarks = filters_dict.get(
                        seq_name, (0, np.zeros(3), (0, 0))
                    )
                    n_hallmarks = n_hallmarks[0]
                    marker_enrichment = f"{marker_enrichment[1]:.4f}"
                    conjugation_genes = conjscan_genes_dict.get(seq_name, [])
                    if len(conjugation_genes):
                        conjugation_genes = ";".join(conjugation_genes)
                    else:
                        conjugation_genes = "NA"
                    amr_genes = amr_genes_dict.get(seq_name, [])
                    if len(amr_genes):
                        amr_genes = ";".join(amr_genes)
                    else:
                        amr_genes = "NA"
                else:
                    marker_enrichment = "NA"
                    n_hallmarks = "NA"
                    conjugation_genes = "NA"
                    amr_genes = "NA"
                fout.write(
                    f"{seq_name}\t{length}\t{terminal_repeat}\t{n_genes}\t{genetic_code}\t{score}\t"
                    f"{fdr}\t{n_hallmarks}\t{marker_enrichment}\t{conjugation_genes}\t{amr_genes}\n"
                )

        # Write virus summary file
        with open(outputs.summary_virus_output, "w") as fout:
            fout.write(
                "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\t"
                "virus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
            )
            for seq_name, score, fdr in itertools.zip_longest(
                virus_name_array, virus_score_array, virus_fdr_array, fillvalue="NA"
            ):
                length = length_dict.get(seq_name, "NA")
                terminal_repeat = terminal_repeat_dict.get(seq_name, "NA")
                coord = provirus_coord_dict.get(seq_name, "NA")
                coord = "-".join(map(str, coord)) if isinstance(coord, tuple) else coord
                n_genes = n_genes_dict.get(seq_name, "NA")
                genetic_code = genetic_code_dict.get(seq_name, "NA")
                score = f"{score:.4f}"
                fdr = f"{fdr:.4f}" if not isinstance(fdr, str) else fdr
                if annotate_exec:
                    _, marker_enrichment, n_hallmarks = filters_dict.get(
                        seq_name, (0, np.zeros(3), (0, 0))
                    )
                    n_hallmarks = n_hallmarks[1]
                    marker_enrichment = f"{marker_enrichment[2]:.4f}"
                    taxonomy = taxonomy_dict.get(seq_name, "Unclassified")
                else:
                    marker_enrichment = "NA"
                    n_hallmarks = "NA"
                    taxonomy = "NA"
                fout.write(
                    f"{seq_name}\t{length}\t{terminal_repeat}\t{coord}\t{n_genes}\t{genetic_code}\t"
                    f"{score}\t{fdr}\t{n_hallmarks}\t{marker_enrichment}\t{taxonomy}\n"
                )

        console.log(
            "Summary files were written to "
            f"[green]{outputs.summary_plasmid_output.name}[/green] and "
            f"[green]{outputs.summary_virus_output.name}[/green]."
        )

        console.log("geNomad summary finished!", style="yellow")
