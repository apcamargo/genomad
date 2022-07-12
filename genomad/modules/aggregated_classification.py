import sys

import numpy as np
from genomad import sequence, utils
from genomad._paths import GenomadOutputs


def branch_attention(w: np.array, b1: np.array, b2: np.array, temperature: float = 2):
    w_1 = np.array(
        [[0.3598502, 2.912244, -1.0668367, 1.3729712, -2.1972055, 0.9363847]]
    )
    w_2 = np.array(
        [[1.5372132, 2.6216774, -2.8225133, 3.0680428, 2.803005, -1.1982375]]
    )
    alpha = np.matmul(w.reshape(-1, 1), w_1) + w_2
    b1 = b1 * alpha[:, 0:3]
    b2 = b2 * alpha[:, 3:6]
    dense_layer_weights = np.array(
        [
            [1.6666023, -1.1003100, -2.1425622],
            [-2.2625937, 2.7540822, -1.5622343],
            [1.9745151, 1.0952991, -2.7467837],
        ]
    )
    dense_layer_bias = np.array([0.14732242, -0.6838019, 0.5594167])
    output = np.matmul((b1 + b2) / 2, dense_layer_weights) + dense_layer_bias
    return utils.softmax(output, temperature)


def main(input_path, output_path, restart, verbose):
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
        output_file=outputs.aggregated_classification_log, verbose=verbose
    )
    # Create a dictionary containing the parameters
    parameter_dict = {}

    # Check if the find-provirus module was executed
    classify_proviruses = utils.check_provirus_execution(
        prefix, input_path, output_path
    )

    # Display the module header
    output_files_path_list = [
        outputs.aggregated_classification_execution_info,
        outputs.aggregated_classification_output,
        outputs.aggregated_classification_npz_output,
    ]
    output_files_description_list = [
        "execution parameters",
        "sequence classification: tabular format",
        "sequence classification: binary format",
    ]
    if classify_proviruses:
        output_files_path_list.extend(
            [
                outputs.provirus_aggregated_classification_output,
                outputs.provirus_aggregated_classification_npz_output,
            ]
        )
        output_files_description_list.extend(
            [
                "provirus classification: tabular format",
                "provirus classification: binary format",
            ]
        )

    utils.display_header(
        console,
        "aggregated-classification",
        (
            "This will aggregate the results of the [cyan]marker-classification[/cyan] "
            "and [cyan]nn-classification[/cyan] modules to classify the input "
            "sequences into chromosome, plasmid, or virus."
        ),
        outputs.aggregated_classification_dir,
        output_files_path_list,
        output_files_description_list,
    )

    # Check if the marker-classification and nn-classification modules were previously executed
    required_inputs = [
        outputs.marker_classification_execution_info,
        outputs.features_npz_output,
        outputs.marker_classification_npz_output,
        outputs.nn_classification_execution_info,
        outputs.nn_classification_npz_output,
    ]
    if classify_proviruses:
        required_inputs.extend(
            [
                outputs.provirus_marker_classification_npz_output,
                outputs.provirus_nn_classification_npz_output,
            ]
        )
    missing_required_inputs = [i.name for i in required_inputs if not i.exists()]
    if len(missing_required_inputs):
        console.error(
            "The following files could not be found: [u]"
            + "[/u], [u]".join(missing_required_inputs)
            + "[/u]. Make sure to execute the [cyan]marker-classification[/cyan] "
            + "and [cyan]nn-classification[/cyan] modules."
        )
        sys.exit(1)

    # Compare the MD5 hashes of the summary, marker-classification, and nn-classification
    # executions
    marker_classification_md5 = utils.get_execution_info(
        outputs.marker_classification_execution_info
    )[0]
    nn_classification_md5 = utils.get_execution_info(
        outputs.nn_classification_execution_info
    )[0]
    input_md5 = utils.get_md5(input_path)
    if (input_md5 != marker_classification_md5) or (input_md5 != nn_classification_md5):
        console.error(
            "Different input FASTA files were used as input for the "
            "[cyan]marker-classification[/cyan], [cyan]nn-classification[/cyan], "
            "and [cyan]aggregated-classification[/cyan] modules. Please execute "
            "all modules using the same input."
        )
        sys.exit(1)

    # Check the input FASTA file
    if not sequence.check_fasta(input_path):
        console.error(
            f"[u]{input_path}[/u] is either empty or contains multiple entries "
            "with the same identifier. Please check your input FASTA file and "
            "execute [cyan]genomad aggregated-classification[/cyan] again."
        )
        sys.exit(1)

    # Print initial log
    console.log("Executing [cyan]genomad aggregated-classification[/cyan].")

    # Check if steps can be skipped
    skip = False
    if (
        outputs.aggregated_classification_execution_info.exists()
        and any(path.exists() for path in output_files_path_list)
        and not restart
    ):
        # Check if the previous execution used the same input file and parameters
        if utils.compare_executions(
            input_path, parameter_dict, outputs.aggregated_classification_execution_info
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
    if not outputs.aggregated_classification_dir.is_dir():
        console.log(
            f"Creating the [green]{outputs.aggregated_classification_dir}[/green] directory."
        )
        outputs.aggregated_classification_dir.mkdir()

    # Write the execution data to `execution_info_output`
    utils.write_execution_info(
        "aggregated_classification",
        input_path,
        parameter_dict,
        outputs.aggregated_classification_execution_info,
    )

    # Compute total marker frequency
    with console.status(
        "Computing the total marker frequencies of the input sequences."
    ):
        contig_marker_freq = np.load(outputs.features_npz_output)["contig_features"]
        contig_marker_freq = contig_marker_freq[:, 15:18].sum(1)
        if classify_proviruses:
            provirus_marker_freq = np.load(outputs.provirus_features_npz_output)[
                "provirus_features"
            ]
            provirus_marker_freq = provirus_marker_freq[:, 15:18].sum(1)
        console.log(
            "The total marker frequencies of the input sequences were computed."
        )

    # Classify sequences and write `aggregated_classification_npz_output`
    if skip and outputs.aggregated_classification_npz_output.exists():
        console.log(
            f"[green]{outputs.aggregated_classification_npz_output.name}[/green] was found. "
            "Skipping sequence classification."
        )
        contig_names = np.load(outputs.aggregated_classification_npz_output)[
            "contig_names"
        ]
        contig_predictions = np.load(outputs.aggregated_classification_npz_output)[
            "predictions"
        ]
    else:
        with console.status("Classifying sequences."):
            contig_names = np.load(outputs.marker_classification_npz_output)[
                "contig_names"
            ]
            marker_predictions = np.load(outputs.marker_classification_npz_output)[
                "predictions"
            ]
            nn_predictions = np.load(outputs.nn_classification_npz_output)[
                "predictions"
            ]
            contig_predictions = branch_attention(
                contig_marker_freq, marker_predictions, nn_predictions
            )
            console.log("Sequences classified.")
        with console.status(
            "Writing sequence classification in binary format to [green]"
            f"{outputs.aggregated_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.aggregated_classification_npz_output,
                contig_names=contig_names,
                predictions=contig_predictions,
            )
            console.log(
                "Sequence classification in binary format written to "
                f"[green]{outputs.aggregated_classification_npz_output.name}[/green]."
            )

    # Write `aggregated_classification_output`
    with console.status(
        "Writing sequence classification in tabular format to [green]"
        f"{outputs.aggregated_classification_output.name}[/green]."
    ):
        with open(outputs.aggregated_classification_output, "w") as fout:
            fout.write("seq_name\tchromosome_score\tplasmid_score\tvirus_score\n")
            for name, scores in zip(contig_names, contig_predictions):
                scores = "".join(map(lambda x: f"{x:.4f}\t", scores)).strip()
                fout.write(f"{name}\t{scores}\n")
        console.log(
            "Sequence classification in tabular format written to [green]"
            f"{outputs.aggregated_classification_output.name}[/green]."
        )

    # Classify proviruses and write `provirus_aggregated_classification_npz_output`
    if (
        skip
        and classify_proviruses
        and outputs.provirus_aggregated_classification_npz_output.exists()
    ):
        console.log(
            f"[green]{outputs.provirus_aggregated_classification_npz_output.name}[/green] was found. "
            "Skipping provirus classification."
        )
        provirus_names = np.load(outputs.provirus_aggregated_classification_npz_output)[
            "provirus_names"
        ]
        provirus_predictions = np.load(
            outputs.provirus_aggregated_classification_npz_output
        )["predictions"]
    elif classify_proviruses:
        with console.status("Classifying proviruses."):
            provirus_names = np.load(outputs.provirus_marker_classification_npz_output)[
                "provirus_names"
            ]
            marker_predictions = np.load(
                outputs.provirus_marker_classification_npz_output
            )["predictions"]
            nn_predictions = np.load(outputs.provirus_nn_classification_npz_output)[
                "predictions"
            ]
            provirus_predictions = branch_attention(
                provirus_marker_freq, marker_predictions, nn_predictions
            )
            console.log("Proviruses classified.")
        with console.status(
            "Writing provirus classification in binary format to [green]"
            f"{outputs.provirus_aggregated_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.provirus_aggregated_classification_npz_output,
                provirus_names=provirus_names,
                predictions=provirus_predictions,
            )
            console.log(
                "Provirus classification in binary format written to "
                f"[green]{outputs.provirus_aggregated_classification_npz_output.name}[/green]."
            )

    # Write `provirus_aggregated_classification_output`
    if classify_proviruses:
        with console.status(
            "Writing provirus classification in tabular format to [green]"
            f"{outputs.provirus_aggregated_classification_output.name}[/green]."
        ):
            with open(outputs.provirus_aggregated_classification_output, "w") as fout:
                fout.write("seq_name\tchromosome_score\tplasmid_score\tvirus_score\n")
                for name, scores in zip(provirus_names, provirus_predictions):
                    scores = "".join(map(lambda x: f"{x:.4f}\t", scores)).strip()
                    fout.write(f"{name}\t{scores}\n")
            console.log(
                "Provirus classification in tabular format written to [green]"
                f"{outputs.provirus_aggregated_classification_output.name}[/green]."
            )

    console.log("geNomad aggregated-classification finished!", style="yellow")
