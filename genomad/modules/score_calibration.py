import sys
from collections import Counter

import numpy as np
from genomad import utils
from genomad._paths import GenomadData, GenomadOutputs


def get_empirical_sample_composition(score_array):
    counts = Counter(score_array.argmax(1))
    counts = [counts[i] for i in range(score_array.shape[1])]
    return np.array(counts) / np.sum(counts)


def score_batch_correction(scores, composition, classifier, weights_file):
    # Apply the smoothing coefficient to shrink the calibration effect in
    # very skewed distributions
    composition_skew = utils.specificity(composition)
    smoothing_coef = 1 - (composition_skew * 0.3)
    uniform_dist = np.ones(3) / 3
    composition = (composition * smoothing_coef) + (uniform_dist * (1 - smoothing_coef))
    # Set the aggregated classifier as the default
    if classifier not in {"marker", "aggregated", "nn"}:
        classifier = "aggregated"
    # Prepare inputs
    composition = np.expand_dims(composition, 0)
    composition = np.repeat(composition, scores.shape[0], 0)
    # Load weights
    npz_file = np.load(weights_file)
    kernel_1 = npz_file[f"kernel_1_{classifier}"]
    bias_1 = npz_file[f"bias_1_{classifier}"]
    kernel_2 = npz_file[f"kernel_2_{classifier}"]
    bias_2 = npz_file[f"bias_2_{classifier}"]
    kernel_3 = npz_file[f"kernel_3_{classifier}"]
    bias_3 = npz_file[f"bias_3_{classifier}"]
    # Inference
    x = np.concatenate([composition, scores], 1)
    x = np.matmul(x, kernel_1) + bias_1
    x = np.tanh(x)
    x = np.matmul(x, kernel_2) + bias_2
    x = np.tanh(x)
    x = np.matmul(x, kernel_3) + bias_3
    return utils.softmax(x)


def write_score_output(output_path, name_array, score_array):
    with open(output_path, "w") as fout:
        fout.write("seq_name\tchromosome_score\tplasmid_score\tvirus_score\n")
        for n, (c_score, p_score, v_score) in zip(name_array, score_array):
            fout.write(f"{n}\t{c_score:.4f}\t{p_score:.4f}\t{v_score:.4f}\n")


def main(input_path, output_path, composition, force_auto, verbose):
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
        output_file=outputs.score_calibration_log, verbose=verbose
    )
    # Create a dictionary containing the parameters
    parameter_dict = {
        "composition": composition,
        "force_auto": force_auto,
    }

    if composition not in {"auto", "metagenome", "virome"}:
        console.error("Invalid value for the [u]composition[/u] parameter.")
        sys.exit(1)

    # Check which modules were previously executed
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

    # Exit if no classification is found
    if not any(
        [
            marker_classification_exec,
            nn_classification_exec,
            aggregated_classification_exec,
        ]
    ):
        console.error(
            "No previous execution of the [cyan]marker-classification[/cyan], "
            "[cyan]nn-classification[/cyan], or [cyan]aggregated-classification[/cyan] "
            "were detected. Please execute at least one of these modules."
        )
        sys.exit(1)

    # Display header
    output_files_path_list = [
        outputs.score_calibration_execution_info,
        outputs.score_calibration_compositions_output,
        outputs.score_calibration_compositions_npz_output,
    ]
    output_files_description_list = [
        "execution parameters",
        "estimated compositions: tabular format",
        "estimated compositions: binary format",
    ]
    if marker_classification_exec:
        output_files_path_list.append(outputs.calibrated_marker_classification_output)
        output_files_path_list.append(
            outputs.calibrated_marker_classification_npz_output
        )
        output_files_description_list.append(
            "calibrated sequence marker-classification data: tabular format"
        )
        output_files_description_list.append(
            "calibrated sequence marker-classification data: binary format"
        )
    if nn_classification_exec:
        output_files_path_list.append(outputs.calibrated_nn_classification_output)
        output_files_path_list.append(outputs.calibrated_nn_classification_npz_output)
        output_files_description_list.append(
            "calibrated sequence nn-classification data: tabular format"
        )
        output_files_description_list.append(
            "calibrated sequence nn-classification data: binary format"
        )
    if aggregated_classification_exec:
        output_files_path_list.append(
            outputs.calibrated_aggregated_classification_output
        )
        output_files_path_list.append(
            outputs.calibrated_aggregated_classification_npz_output
        )
        output_files_description_list.append(
            "calibrated sequence aggregated-classification data: tabular format"
        )
        output_files_description_list.append(
            "calibrated sequence aggregated-classification data: binary format"
        )
    if marker_classification_provirus_exec:
        output_files_path_list.append(
            outputs.provirus_calibrated_marker_classification_output
        )
        output_files_path_list.append(
            outputs.provirus_calibrated_marker_classification_npz_output
        )
        output_files_description_list.append(
            "calibrated provirus marker-classification data: tabular format"
        )
        output_files_description_list.append(
            "calibrated provirus marker-classification data: binary format"
        )
    if nn_classification_provirus_exec:
        output_files_path_list.append(
            outputs.provirus_calibrated_nn_classification_output
        )
        output_files_path_list.append(
            outputs.provirus_calibrated_nn_classification_npz_output
        )
        output_files_description_list.append(
            "calibrated provirus nn-classification data: tabular format"
        )
        output_files_description_list.append(
            "calibrated provirus nn-classification data: binary format"
        )
    if aggregated_classification_provirus_exec:
        output_files_path_list.append(
            outputs.provirus_calibrated_aggregated_classification_output
        )
        output_files_path_list.append(
            outputs.provirus_calibrated_aggregated_classification_npz_output
        )
        output_files_description_list.append(
            "calibrated provirus aggregated-classification data: tabular format"
        )
        output_files_description_list.append(
            "calibrated provirus aggregated-classification data: binary format"
        )
    utils.display_header(
        console,
        "score-calibration",
        (
            "This will calibrate the classification scores based on the sample composition."
        ),
        outputs.score_calibration_dir,
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
    if len(set(md5_list)) > 1:
        console.error(
            "Different input FASTA files were used as input for the different modules. Please "
            "execute all the modules using the same input."
        )
        sys.exit(1)

    # Print initial log
    console.log("Executing [cyan]genomad score-calibration[/cyan].")

    # If necessary, create the subdirectory
    if not outputs.score_calibration_dir.is_dir():
        console.log(
            f"Creating the [green]{outputs.score_calibration_dir}[/green] directory."
        )
        outputs.score_calibration_dir.mkdir()

    # Write the execution data to `execution_info_output`
    utils.write_execution_info(
        "score_calibration",
        input_path,
        parameter_dict,
        outputs.score_calibration_execution_info,
    )

    with console.status("Loading the uncalibrated classification scores."):
        # Load the scores and names into the `score_dict` dictionary
        score_dict = {}
        if marker_classification_exec:
            score_dict["contig_names"] = np.load(
                outputs.marker_classification_npz_output
            )["contig_names"]
            score_dict["marker_scores"] = np.load(
                outputs.marker_classification_npz_output
            )["predictions"]
        if nn_classification_exec:
            score_dict["contig_names"] = np.load(outputs.nn_classification_npz_output)[
                "contig_names"
            ]
            score_dict["nn_scores"] = np.load(outputs.nn_classification_npz_output)[
                "predictions"
            ]
        if aggregated_classification_exec:
            score_dict["contig_names"] = np.load(
                outputs.aggregated_classification_npz_output
            )["contig_names"]
            score_dict["aggregated_scores"] = np.load(
                outputs.aggregated_classification_npz_output
            )["predictions"]
        if marker_classification_provirus_exec:
            score_dict["provirus_names"] = np.load(
                outputs.provirus_marker_classification_npz_output
            )["provirus_names"]
            score_dict["provirus_marker_scores"] = np.load(
                outputs.provirus_marker_classification_npz_output
            )["predictions"]
        if nn_classification_provirus_exec:
            score_dict["provirus_names"] = np.load(
                outputs.provirus_nn_classification_npz_output
            )["provirus_names"]
            score_dict["provirus_nn_scores"] = np.load(
                outputs.provirus_nn_classification_npz_output
            )["predictions"]
        if aggregated_classification_provirus_exec:
            score_dict["provirus_names"] = np.load(
                outputs.provirus_aggregated_classification_npz_output
            )["provirus_names"]
            score_dict["provirus_aggregated_scores"] = np.load(
                outputs.provirus_aggregated_classification_npz_output
            )["predictions"]
        console.log("The uncalibrated classification scores were loaded.")

    # Get sample composition
    n_sequences = (
        len(score_dict["contig_names"]) + len(score_dict["provirus_names"])
        if find_proviruses_exec
        else len(score_dict["contig_names"])
    )
    if n_sequences < 1_000 and composition == "auto" and not force_auto:
        console.warning(
            "Your sample has less than 1,000 sequences, which does not allow precise "
            "composition estimation. The 'metagenome' preset will be used instead. If "
            "you want to force the empirical sample composition estimation anyway, "
            "use the [cyan]--force-auto[/cyan] option.",
        )
        composition = "metagenome"

    # Get the empirical sample composition
    if composition == "auto":
        composition_dict = {}
        if marker_classification_exec and marker_classification_provirus_exec:
            composition_dict["marker"] = get_empirical_sample_composition(
                np.concatenate(
                    [score_dict["marker_scores"], score_dict["provirus_marker_scores"]]
                )
            )
        elif marker_classification_exec:
            composition_dict["marker"] = get_empirical_sample_composition(
                score_dict["marker_scores"]
            )
        if nn_classification_exec and nn_classification_provirus_exec:
            composition_dict["nn"] = get_empirical_sample_composition(
                np.concatenate(
                    [score_dict["nn_scores"], score_dict["provirus_nn_scores"]]
                )
            )
        elif nn_classification_exec:
            composition_dict["nn"] = get_empirical_sample_composition(
                score_dict["nn_scores"]
            )
        if aggregated_classification_exec and aggregated_classification_provirus_exec:
            composition_dict["aggregated"] = get_empirical_sample_composition(
                np.concatenate(
                    [
                        score_dict["aggregated_scores"],
                        score_dict["provirus_aggregated_scores"],
                    ]
                )
            )
        elif aggregated_classification_exec:
            composition_dict["aggregated"] = get_empirical_sample_composition(
                score_dict["aggregated_scores"]
            )
    # Use presets
    elif composition == "metagenome":
        composition_dict = {
            "marker": np.array([0.84, 0.05, 0.11]),
            "nn": np.array([0.67, 0.20, 0.13]),
            "aggregated": np.array([0.72, 0.17, 0.11]),
        }
    elif composition == "virome":
        composition_dict = {
            "marker": np.array([0.26, 0.004, 0.736]),
            "nn": np.array([0.23, 0.06, 0.71]),
            "aggregated": np.array([0.24, 0.025, 0.735]),
        }

    with console.status(
        "Writing estimated compositions to [green]"
        f"{outputs.score_calibration_compositions_output.name}[/green] and [green]"
        f"{outputs.score_calibration_compositions_npz_output.name}[/green]."
    ):
        np.savez_compressed(
            outputs.score_calibration_compositions_npz_output,
            marker=composition_dict.get("marker", np.array([0.0, 0.0, 0.0])),
            nn=composition_dict.get("nn", np.array([0.0, 0.0, 0.0])),
            aggregated=composition_dict.get("aggregated", np.array([0.0, 0.0, 0.0])),
        )
        with open(outputs.score_calibration_compositions_output, "w") as fout:
            fout.write("model\tchromosome\tplasmid\tvirus\n")
            for k, v in composition_dict.items():
                v = "\t".join([f"{i:.4f}" for i in v])
                fout.write(f"{k}\t{v}\n")
        console.log(
            "Estimated compositions written to [green]"
            f"{outputs.score_calibration_compositions_output.name}[/green] and [green]"
            f"{outputs.score_calibration_compositions_npz_output.name}[/green]."
        )

    # Perform score calibration
    with console.status("Calibrating scores."):
        calibrated_score_dict = {}
        for k, v in score_dict.items():
            if "marker_scores" in k:
                calibrated_score_dict[k] = score_batch_correction(
                    v,
                    composition_dict["marker"],
                    "marker",
                    GenomadData.score_calibration_weights_file,
                )
            elif "nn_scores" in k:
                calibrated_score_dict[k] = score_batch_correction(
                    v,
                    composition_dict["nn"],
                    "nn",
                    GenomadData.score_calibration_weights_file,
                )
            elif "aggregated_scores" in k:
                calibrated_score_dict[k] = score_batch_correction(
                    v,
                    composition_dict["aggregated"],
                    "aggregated",
                    GenomadData.score_calibration_weights_file,
                )
        console.log("The calibrated scores were computed.")

    # Write the calibrated scores to output files
    if marker_classification_exec:
        with console.status(
            "Writing calibrated scores in binary format to [green]"
            f"{outputs.calibrated_marker_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.calibrated_marker_classification_npz_output,
                contig_names=score_dict["contig_names"],
                predictions=calibrated_score_dict["marker_scores"],
            )
            console.log(
                "Calibrated scores in binary format written to "
                f"[green]{outputs.calibrated_marker_classification_npz_output.name}[/green]."
            )
        with console.status(
            "Writing calibrated scores in tabular format to [green]"
            f"{outputs.calibrated_marker_classification_output.name}[/green]."
        ):
            write_score_output(
                outputs.calibrated_marker_classification_output,
                score_dict["contig_names"],
                calibrated_score_dict["marker_scores"],
            )
            console.log(
                "Calibrated scores in tabular format written to "
                f"[green]{outputs.calibrated_marker_classification_output.name}[/green]."
            )
    if marker_classification_provirus_exec:
        with console.status(
            "Writing calibrated scores in binary format to [green]"
            f"{outputs.provirus_calibrated_marker_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.provirus_calibrated_marker_classification_npz_output,
                provirus_names=score_dict["provirus_names"],
                predictions=calibrated_score_dict["provirus_marker_scores"],
            )
            console.log(
                "Calibrated scores in binary format written to "
                f"[green]{outputs.provirus_calibrated_marker_classification_npz_output.name}[/green]."
            )
        with console.status(
            "Writing calibrated scores in tabular format to [green]"
            f"{outputs.provirus_calibrated_marker_classification_output.name}[/green]."
        ):
            write_score_output(
                outputs.provirus_calibrated_marker_classification_output,
                score_dict["provirus_names"],
                calibrated_score_dict["provirus_marker_scores"],
            )
            console.log(
                "Calibrated scores in tabular format written to "
                f"[green]{outputs.provirus_calibrated_marker_classification_output.name}[/green]."
            )
    if nn_classification_exec:
        with console.status(
            "Writing calibrated scores in binary format to [green]"
            f"{outputs.calibrated_nn_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.calibrated_nn_classification_npz_output,
                contig_names=score_dict["contig_names"],
                predictions=calibrated_score_dict["nn_scores"],
            )
            console.log(
                "Calibrated scores in binary format written to "
                f"[green]{outputs.calibrated_nn_classification_npz_output.name}[/green]."
            )
        with console.status(
            "Writing calibrated scores in tabular format to [green]"
            f"{outputs.calibrated_nn_classification_output.name}[/green]."
        ):
            write_score_output(
                outputs.calibrated_nn_classification_output,
                score_dict["contig_names"],
                calibrated_score_dict["nn_scores"],
            )
            console.log(
                "Calibrated scores in tabular format written to "
                f"[green]{outputs.calibrated_nn_classification_output.name}[/green]."
            )
    if nn_classification_provirus_exec:
        with console.status(
            "Writing calibrated scores in binary format to [green]"
            f"{outputs.provirus_calibrated_nn_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.provirus_calibrated_nn_classification_npz_output,
                provirus_names=score_dict["provirus_names"],
                predictions=calibrated_score_dict["provirus_nn_scores"],
            )
            console.log(
                "Calibrated scores in binary format written to "
                f"[green]{outputs.provirus_calibrated_nn_classification_npz_output.name}[/green]."
            )
        with console.status(
            "Writing calibrated scores in tabular format to [green]"
            f"{outputs.provirus_calibrated_nn_classification_output.name}[/green]."
        ):
            write_score_output(
                outputs.provirus_calibrated_nn_classification_output,
                score_dict["provirus_names"],
                calibrated_score_dict["provirus_nn_scores"],
            )
            console.log(
                "Calibrated scores in tabular format written to "
                f"[green]{outputs.provirus_calibrated_nn_classification_output.name}[/green]."
            )
    if aggregated_classification_exec:
        with console.status(
            "Writing calibrated scores in binary format to [green]"
            f"{outputs.calibrated_aggregated_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.calibrated_aggregated_classification_npz_output,
                contig_names=score_dict["contig_names"],
                predictions=calibrated_score_dict["aggregated_scores"],
            )
            console.log(
                "Calibrated scores in binary format written to "
                f"[green]{outputs.calibrated_aggregated_classification_npz_output.name}[/green]."
            )
        with console.status(
            "Writing calibrated scores in tabular format to [green]"
            f"{outputs.calibrated_aggregated_classification_output.name}[/green]."
        ):
            write_score_output(
                outputs.calibrated_aggregated_classification_output,
                score_dict["contig_names"],
                calibrated_score_dict["aggregated_scores"],
            )
            console.log(
                "Calibrated scores in tabular format written to "
                f"[green]{outputs.calibrated_aggregated_classification_output.name}[/green]."
            )
    if aggregated_classification_provirus_exec:
        with console.status(
            "Writing calibrated scores in binary format to [green]"
            f"{outputs.provirus_calibrated_aggregated_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.provirus_calibrated_aggregated_classification_npz_output,
                provirus_names=score_dict["provirus_names"],
                predictions=calibrated_score_dict["provirus_aggregated_scores"],
            )
            console.log(
                "Calibrated scores in binary format written to "
                f"[green]{outputs.provirus_calibrated_aggregated_classification_npz_output.name}[/green]."
            )
        with console.status(
            "Writing calibrated scores in tabular format to [green]"
            f"{outputs.provirus_calibrated_aggregated_classification_output.name}[/green]."
        ):
            write_score_output(
                outputs.provirus_calibrated_aggregated_classification_output,
                score_dict["provirus_names"],
                calibrated_score_dict["provirus_aggregated_scores"],
            )
            console.log(
                "Calibrated scores in tabular format written to "
                f"[green]{outputs.provirus_calibrated_aggregated_classification_output.name}[/green]."
            )
    console.log("geNomad score-calibration finished!", style="yellow")
