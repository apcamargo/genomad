import shutil
import sys

from genomad import database, mmseqs2, prodigal, sequence, taxonomy, utils
from genomad._paths import GenomadOutputs


def write_genes_output(genes_output, database_obj, prodigal_obj, mmseqs2_obj):
    marker_annotation_dict = database_obj.get_marker_annotation()
    gene_match_dict = mmseqs2_obj.get_matches()
    taxdb = database_obj.get_taxdb()
    with open(genes_output, "w") as fout:
        fout.write(
            "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\tmarker\t"
            "evalue\tbitscore\tuscg\tplasmid_hallmark\tvirus_hallmark\ttaxid\ttaxname\t"
            "annotation_conjscan\tannotation_amr\tannotation_accessions\tannotation_description\n"
        )
        for (
            contig,
            gene_num,
            start,
            end,
            strand,
            rbs,
            code,
            gc,
        ) in prodigal_obj.proteins():
            gene = f"{contig}_{gene_num}"
            match, evalue, bitscore, taxid = gene_match_dict.get(
                gene, ("NA", "NA", "NA", 1)
            )
            taxname = taxdb.taxid2name[taxid] if taxid != 1 else "NA"
            (
                uscg,
                plasmid_hallmark,
                virus_hallmark,
                conjscan,
                amr,
                accession,
                description,
            ) = marker_annotation_dict.get(match, (0, 0, 0, "NA", "NA", "NA", "NA"))
            gene_length = end - start + 1
            fout.write(
                f"{gene}\t{start}\t{end}\t{gene_length}\t{strand}\t{gc:.3f}\t{code}\t{rbs}\t"
                f"{match}\t{evalue}\t{bitscore}\t{uscg}\t{plasmid_hallmark}\t{virus_hallmark}\t"
                f"{taxid}\t{taxname}\t{conjscan}\t{amr}\t{accession}\t{description}\n"
            )


def main(
    input_path,
    output_path,
    database_path,
    use_minimal_db,
    restart,
    threads,
    verbose,
    sensitivity,
    evalue,
    splits,
    cleanup,
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
    console = utils.HybridConsole(output_file=outputs.annotate_log, verbose=verbose)
    # Create a dictionary containing the parameters
    parameter_dict = {
        "use_minimal_db": use_minimal_db,
        "sensitivity": sensitivity,
        "evalue": evalue,
    }

    # Display the module header
    utils.display_header(
        console,
        "annotate",
        (
            "This will perform gene calling in the input sequences and "
            "annotate the predicted proteins with geNomad's markers."
        ),
        outputs.annotate_dir,
        [
            outputs.annotate_execution_info,
            outputs.annotate_genes_output,
            outputs.annotate_taxonomy_output,
            outputs.annotate_mmseqs2_output,
            outputs.annotate_proteins_output,
        ],
        [
            "execution parameters",
            "gene annotation data",
            "taxonomic assignment",
            "MMseqs2 output file",
            "protein FASTA file",
        ],
    )

    # Check if the required binaries are in the PATH
    required_executables = ["mmseqs", "prodigal-gv"]
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
            "execute [cyan]genomad annotate[/cyan] again."
        )
        sys.exit(1)

    # Print initial log
    console.log("Executing [cyan]genomad annotate[/cyan].")

    # Check if steps can be skipped
    skip = False
    if (
        outputs.annotate_execution_info.exists()
        and any(
            [
                outputs.annotate_proteins_output.exists(),
                outputs.annotate_genes_output.exists(),
            ]
        )
        and not restart
    ):
        # Check if the previous execution used the same input file and parameters
        if utils.compare_executions(
            input_path, parameter_dict, outputs.annotate_execution_info
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
    if not outputs.annotate_dir.is_dir():
        console.log(f"Creating the [green]{outputs.annotate_dir}[/green] directory.")
        outputs.annotate_dir.mkdir()

    # Write the execution data to `execution_info_output`
    utils.write_execution_info(
        "annotate", input_path, parameter_dict, outputs.annotate_execution_info
    )

    # Initialize the database object
    database_obj = database.Database(database_path)

    # Run prodigal-gv
    prodigal_obj = prodigal.Prodigal(input_path, outputs.annotate_proteins_output)
    if skip and outputs.annotate_proteins_output.exists():
        console.log(
            f"[green]{outputs.annotate_proteins_output.name}[/green] was found. "
            "Skipping gene prediction with prodigal-gv."
        )
    else:
        with console.status("Predicting proteins with prodigal-gv."):
            prodigal_obj.run_parallel_prodigal(threads)
            console.log(
                "Proteins predicted with prodigal-gv were written to "
                f"[green]{outputs.annotate_proteins_output.name}[/green]."
            )

    # Run MMseqs2
    mmseqs2_obj = mmseqs2.MMseqs2(
        outputs.annotate_mmseqs2_output,
        outputs.annotate_mmseqs2_dir,
        outputs.annotate_proteins_output,
        database_obj,
        use_minimal_db=use_minimal_db,
    )
    if skip and outputs.annotate_mmseqs2_output.exists():
        console.log(
            f"[green]{outputs.annotate_mmseqs2_output.name}[/green] was found. "
            "Skipping protein annotation with MMseqs2."
        )
    else:
        with console.status(
            "Annotating proteins with MMseqs2 and geNomad database "
            f"(v{database_obj.version})."
        ):
            mmseqs2_obj.run_mmseqs2(threads, sensitivity, evalue, splits)
            console.log(
                "Proteins annotated with MMseqs2 and geNomad database "
                f"(v{database_obj.version}) were written to "
                f"[green]{outputs.annotate_mmseqs2_output.name}[/green]."
            )
    if cleanup and outputs.annotate_mmseqs2_dir.is_dir():
        console.log(f"Deleting [green]{outputs.annotate_mmseqs2_dir.name}[/green].")
        shutil.rmtree(outputs.annotate_mmseqs2_dir)

    # Write `annotate_genes_output`
    with console.status(
        f"Writing [green]{outputs.annotate_genes_output.name}[/green]."
    ):
        write_genes_output(
            outputs.annotate_genes_output, database_obj, prodigal_obj, mmseqs2_obj
        )
        console.log(
            f"Gene data was written to [green]{outputs.annotate_genes_output.name}[/green]."
        )

    # Write `annotate_taxonomy_output`
    with console.status(
        f"Writing [green]{outputs.annotate_taxonomy_output.name}[/green]."
    ):
        taxonomy.write_taxonomic_assignment(
            outputs.annotate_taxonomy_output,
            outputs.annotate_genes_output,
            database_obj,
        )
        console.log(
            f"Taxonomic assignment data was written to [green]{outputs.annotate_taxonomy_output.name}[/green]."
        )

    console.log("geNomad annotate finished!", style="yellow")
