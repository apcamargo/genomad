import os
import shutil
import sys
from pathlib import Path

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import numpy as np
from genomad import sequence, utils
from genomad._paths import GenomadData, GenomadOutputs


def main(input_path, output_path, single_window, batch_size, restart, verbose, cleanup):
    # To avoid having other modules lagging due to the slow TensorFlow import,
    # the `tensorflow` and `genomad.neural_network` modules are loaded inside `main`.
    # Additionally, the following functions that use the `tensorflow` module are
    # defined within `main`: `write_tfrecord`, `generate_data`, `prepare_sample`,
    # `parse_tfrecord`, `get_dataset`
    import tensorflow as tf
    from genomad import neural_network

    AUTOTUNE = tf.data.experimental.AUTOTUNE

    def write_tfrecord(data: np.array, tfrecords_dir: Path, window_count: int):
        tfrecord_path = tfrecords_dir / f"{window_count}.tfrec"
        tfrecord_path = str(tfrecord_path)
        with tf.io.TFRecordWriter(tfrecord_path) as writer:
            for row in data:
                feature = tf.train.Feature(int64_list=tf.train.Int64List(value=row))
                example = tf.train.Example(
                    features=tf.train.Features(feature={"sequence": feature})
                )
                writer.write(example.SerializeToString())

    def generate_data(
        fasta_path: Path,
        tfrecords_dir: Path,
        single_window: bool = False,
        n_records_per_file: int = 10_000,
    ):
        contig_id_array = []
        contig_name_array = []
        tokenized_seq_data = []
        window_count = 0
        max_windows = 1 if single_window else None
        for contig_id, seq in enumerate(sequence.read_fasta(fasta_path, strip_n=True)):
            contig_name_array.append(seq.accession)
            for seq_window in sequence.seq_windows(
                seq, 6_000, 2_500, max_windows=max_windows
            ):
                tokenized_seq_window = seq_window.seq_ascii.ljust(6_000, b"N")
                tokenized_seq_window = sequence.tokenize_dna(tokenized_seq_window, 4)
                tokenized_seq_data.append(tokenized_seq_window)
                contig_id_array.append(contig_id)
                window_count += 1
                if window_count % n_records_per_file == 0:
                    write_tfrecord(tokenized_seq_data, tfrecords_dir, window_count)
                    tokenized_seq_data = []
        if len(tokenized_seq_data):
            write_tfrecord(tokenized_seq_data, tfrecords_dir, window_count)
        return np.array(contig_name_array), np.array(contig_id_array)

    def prepare_sample(features):
        return features["sequence"]

    def parse_tfrecord(example):
        feature_description = {
            "sequence": tf.io.FixedLenFeature([5_997], tf.int64),
        }
        return tf.io.parse_single_example(example, feature_description)

    def get_dataset(filenames, batch_size):
        return (
            tf.data.TFRecordDataset(filenames, num_parallel_reads=None)
            .map(parse_tfrecord, num_parallel_calls=AUTOTUNE, deterministic=True)
            .map(prepare_sample, num_parallel_calls=AUTOTUNE, deterministic=True)
            .batch(batch_size)
            .prefetch(AUTOTUNE)
        )

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
        output_file=outputs.nn_classification_log, verbose=verbose
    )
    # Create a dictionary containing the parameters
    parameter_dict = {
        "single_window": single_window,
    }
    # Check if the find-provirus module was executed
    classify_proviruses = utils.check_provirus_execution(
        prefix, input_path, output_path
    )

    # Display the module header
    output_files_path_list = [
        outputs.nn_classification_execution_info,
        outputs.encoded_sequences_dir,
        outputs.nn_classification_output,
        outputs.nn_classification_npz_output,
    ]
    output_files_description_list = [
        "execution parameters",
        "directory containing encoded sequence data",
        "contig classification: tabular format",
        "contig classification: binary format",
    ]
    if classify_proviruses:
        output_files_path_list.extend(
            [
                outputs.encoded_proviruses_dir,
                outputs.provirus_nn_classification_output,
                outputs.provirus_nn_classification_npz_output,
            ]
        )
        output_files_description_list.extend(
            [
                "directory containing encoded sequence data",
                "provirus classification: tabular format",
                "provirus classification: binary format",
            ]
        )
    utils.display_header(
        console,
        "nn-classification",
        (
            "This will classify the input sequences into chromosome, plasmid, "
            "or virus based on the nucleotide sequence."
        ),
        outputs.nn_classification_dir,
        output_files_path_list,
        output_files_description_list,
    )

    # Check the input FASTA file
    if not sequence.check_fasta(input_path):
        console.error(
            f"[u]{input_path}[/u] is either empty or contains multiple entries "
            "with the same identifier. Please check your input FASTA file and "
            "execute [cyan]genomad nn-classification[/cyan] again."
        )
        sys.exit(1)

    # Print initial log
    console.log("Executing [cyan]genomad nn-classification[/cyan].")

    # Check if steps can be skipped
    skip = False
    if (
        outputs.nn_classification_execution_info.exists()
        and any(path.exists() for path in output_files_path_list)
        and not restart
    ):
        # Check if the previous execution used the same input file and parameters
        if utils.compare_executions(
            input_path, parameter_dict, outputs.nn_classification_execution_info
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
    if not outputs.nn_classification_dir.is_dir():
        console.log(
            f"Creating the [green]{outputs.nn_classification_dir}[/green] directory."
        )
        outputs.nn_classification_dir.mkdir()

    # Write the execution data to `execution_info_output`
    utils.write_execution_info(
        "nn_classification",
        input_path,
        parameter_dict,
        outputs.nn_classification_execution_info,
    )

    # Encode sequences, write `.tfrec` files, and write `seq_window_id_output`
    if (
        skip
        and outputs.seq_window_id_output.exists()
        and len(list(outputs.encoded_sequences_dir.glob("*.tfrec")))
    ):
        console.log(
            f"[green]{outputs.encoded_sequences_dir.name}[/green] was found. "
            "Skipping sequence encoding."
        )
        contig_names = np.load(outputs.seq_window_id_output)["contig_names"]
        contig_ids = np.load(outputs.seq_window_id_output)["contig_ids"]
    else:
        if outputs.encoded_sequences_dir.is_dir():
            shutil.rmtree(outputs.encoded_sequences_dir)
        console.log(
            f"Creating the [green]{outputs.encoded_sequences_dir}[/green] directory."
        )
        outputs.encoded_sequences_dir.mkdir()
        with console.status("Encoding sequences."):
            contig_names, contig_ids = generate_data(
                input_path, outputs.encoded_sequences_dir, single_window
            )
            np.savez_compressed(
                outputs.seq_window_id_output,
                contig_names=contig_names,
                contig_ids=contig_ids,
            )
            console.log(
                "Encoded sequence data written to "
                f"[green]{outputs.encoded_sequences_dir.name}[/green]."
            )

    # Encode proviruses, write `.tfrec` files, and write `provirus_window_id_output`
    if (
        skip
        and classify_proviruses
        and outputs.provirus_window_id_output.exists()
        and len(list(outputs.encoded_proviruses_dir.glob("*.tfrec")))
    ):
        console.log(
            f"[green]{outputs.encoded_proviruses_dir.name}[/green] was found. "
            "Skipping provirus encoding."
        )
        provirus_names = np.load(outputs.provirus_window_id_output)["provirus_names"]
        provirus_ids = np.load(outputs.provirus_window_id_output)["provirus_ids"]
    elif classify_proviruses:
        if outputs.encoded_proviruses_dir.is_dir():
            shutil.rmtree(outputs.encoded_proviruses_dir)
        console.log(
            f"Creating the [green]{outputs.encoded_proviruses_dir}[/green] directory."
        )
        outputs.encoded_proviruses_dir.mkdir()
        with console.status("Encoding proviruses."):
            provirus_names, provirus_ids = generate_data(
                outputs.find_proviruses_nucleotide_output,
                outputs.encoded_proviruses_dir,
                single_window,
            )
            np.savez_compressed(
                outputs.provirus_window_id_output,
                provirus_names=provirus_names,
                provirus_ids=provirus_ids,
            )
            console.log(
                "Encoded provirus data written to "
                f"[green]{outputs.encoded_proviruses_dir.name}[/green]."
            )

    # Classify sequences and write `nn_classification_npz_output`
    if skip and outputs.nn_classification_npz_output.exists():
        console.log(
            f"[green]{outputs.nn_classification_npz_output.name}[/green] was found. Skipping "
            "sequence classification."
        )
        contig_names = np.load(outputs.nn_classification_npz_output)["contig_names"]
        contig_predictions = np.load(outputs.nn_classification_npz_output)[
            "predictions"
        ]
    else:
        tfrecord_list = utils.natsort(
            tf.io.gfile.glob(f"{outputs.encoded_sequences_dir}/*.tfrec")
        )
        if not len(tfrecord_list):
            console.error("No sequences were found. Please check your input FASTA.")
            sys.exit(1)
        with console.status("Classifying sequences."):
            nn_model = neural_network.create_classifier()
            nn_model.load_weights(GenomadData.nn_model_file)
            contig_predictions = nn_model.predict(
                get_dataset(tfrecord_list, batch_size), verbose=0
            )
            contig_predictions = tf.math.segment_mean(contig_predictions, contig_ids)
            console.log("Sequences classified.")
        with console.status(
            "Writing sequence classification in binary format to [green]"
            f"{outputs.nn_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.nn_classification_npz_output,
                contig_names=contig_names,
                predictions=contig_predictions,
            )
            console.log(
                "Sequence classification in binary format written to [green]"
                f"{outputs.nn_classification_npz_output.name}[/green]."
            )
    if cleanup and outputs.encoded_sequences_dir.is_dir():
        console.log("Deleting encoded sequence data.")
        shutil.rmtree(outputs.encoded_sequences_dir)

    # Write `nn_classification_output`
    with console.status(
        "Writing sequence classification in tabular format to [green]"
        f"{outputs.nn_classification_output.name}[/green]."
    ):
        with open(outputs.nn_classification_output, "w") as fout:
            fout.write("seq_name\tchromosome_score\tplasmid_score\tvirus_score\n")
            for name, scores in zip(contig_names, contig_predictions):
                scores = "".join(map(lambda x: f"{x:.4f}\t", scores)).strip()
                fout.write(f"{name}\t{scores}\n")
        console.log(
            "Sequence classification in tabular format written to [green]"
            f"{outputs.nn_classification_output.name}[/green]."
        )

    # Classify proviruses and write `provirus_nn_classification_npz_output`
    if skip and outputs.provirus_nn_classification_npz_output.exists():
        console.log(
            f"[green]{outputs.provirus_nn_classification_npz_output.name}[/green] was found. "
            "Skipping provirus classification."
        )
        provirus_names = np.load(outputs.provirus_nn_classification_npz_output)[
            "provirus_names"
        ]
        provirus_predictions = np.load(outputs.provirus_nn_classification_npz_output)[
            "predictions"
        ]
    elif classify_proviruses:
        with console.status("Classifying proviruses."):
            nn_model = neural_network.create_classifier()
            nn_model.load_weights(GenomadData.nn_model_file)
            tfrecord_list = utils.natsort(
                tf.io.gfile.glob(f"{outputs.encoded_proviruses_dir}/*.tfrec")
            )
            provirus_predictions = nn_model.predict(
                get_dataset(tfrecord_list, batch_size), verbose=0
            )
            provirus_predictions = tf.math.segment_mean(
                provirus_predictions, provirus_ids
            )
            console.log("Proviruses classified.")
        with console.status(
            "Writing provirus classification in binary format to [green]"
            f"{outputs.provirus_nn_classification_npz_output.name}[/green]."
        ):
            np.savez_compressed(
                outputs.provirus_nn_classification_npz_output,
                provirus_names=provirus_names,
                predictions=provirus_predictions,
            )
            console.log(
                "Provirus classification in binary format written to [green]"
                f"{outputs.provirus_nn_classification_npz_output.name}[/green]."
            )
    if cleanup and outputs.encoded_proviruses_dir.is_dir():
        console.log("Deleting encoded provirus data.")
        shutil.rmtree(outputs.encoded_proviruses_dir)

    # Write `provirus_nn_classification_output`
    if classify_proviruses:
        with console.status(
            "Writing provirus classification in tabular format to [green]"
            f"{outputs.provirus_nn_classification_output.name}[/green]."
        ):
            with open(outputs.provirus_nn_classification_output, "w") as fout:
                fout.write("seq_name\tchromosome_score\tplasmid_score\tvirus_score\n")
                for name, scores in zip(provirus_names, provirus_predictions):
                    scores = "".join(map(lambda x: f"{x:.4f}\t", scores)).strip()
                    fout.write(f"{name}\t{scores}\n")
            console.log(
                "Provirus classification in tabular format written to [green]"
                f"{outputs.provirus_nn_classification_output.name}[/green]."
            )

    console.log("geNomad nn-classification finished!", style="yellow")
