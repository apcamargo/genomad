# geNomad

geNomad: Identification of mobile genetic elements

## Features

geNomad's primary goal is to identify viruses and plasmids in sequencing data (isolates, metagenomes, and metatranscriptomes). It also provides a couple of additional features that can help you in your analysis:

- Taxonomic assignment of viral genomes.
- Identification of viruses integrated in host genomes (proviruses).
- Functional annotation of proteins.

## Documentation

For installation instructions, information about how geNomad works, and a detailed explanation of how to execute it, please check the full documentation: https://portal.nersc.gov/genomad/

## Quick start

We recommend users to read the [documentation](https://portal.nersc.gov/genomad/) before starting to use geNomad. If you are in a rush, however, you can follow this quick step-by-step example.

### Installation

First, you need to install geNomad. There's a couple of ways to do that, but here we will use [conda](https://docs.conda.io/en/latest/) as it will handle all dependencies for us.

```
# Create a conda environment for geNomad
conda create -n genomad -c conda-forge -c bioconda genomad
# Activate the geNomad environment
conda activate genomad
```

### Downloading the database

geNomad depends on a database that contains the profiles of the markers that are used to classify sequences, their taxonomic information, their functional annotation, etc. So, you should first download the database to your current directory:

```
genomad download-database .
```

The database will be contained within the `genomad_db` directory.

### Executing geNomad

Now you are ready to go! geNomad works by executing a series of modules sequentially (more on that in the documentation), but we provide a convenient `end-to-end` command that will execute the entire pipeline for you in one go.

In this example, we will use an *Escherichia coli* genome ([GCF_000008865.2](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000008865.2/)) as input. You can use any FASTA file containing nucleotide sequences as input. geNomad will work for isolate genomes, metagenomes, and metatranscriptomes.

The command to execute geNomad is structured like this:

```
genomad end-to-end [OPTIONS] INPUT OUTPUT DATABASE
```

So, to run the full geNomad pipeline (`end-to-end` command), taking a nucleotide FASTA file (`GCF_000008865.2.fna.gz`) and the database (`genomad_db`) as input, we will execute the following command:

```
genomad end-to-end --cleanup --splits 16 GCF_000008865.2.fna.gz genomad_output genomad_db
```

The results will be written inside the `genomad_output` directory.

Three important details about the command above:

- The `--cleanup` option was used to force geNomad to delete intermediate files that were generated during the execution. This will save you some storage space.
- The `--splits 16` parameter was used here to make it possible to run this example in a notebook. geNomad searches a big database of protein profiles that take up a lot of space in memory. To prevent the execution from failing due to insufficient memory, we can use the `--splits` parameter to split the seach into chuncks. If you are running geNomad in a big server you might not need to split your search, increasing the execution speed.
- Note that the input FASTA file that I used as input was compressed. This is possible because geNomad supports input files compressed as `.gz`, `.bz2`, or `.xz`.

### Understanding the outputs

In this example, the results of geNomad's analysis will be written to the `genomad_output` directory, which will look like this:

```
genomad_output
├── GCF_000008865.2_aggregated_classification
├── GCF_000008865.2_aggregated_classification.log
├── GCF_000008865.2_annotate
├── GCF_000008865.2_annotate.log
├── GCF_000008865.2_find_proviruses
├── GCF_000008865.2_find_proviruses.log
├── GCF_000008865.2_marker_classification
├── GCF_000008865.2_marker_classification.log
├── GCF_000008865.2_nn_classification
├── GCF_000008865.2_nn_classification.log
├── GCF_000008865.2_summary
╰── GCF_000008865.2_summary.log
```

As mentioned above, geNomad works by executing several modules sequentially. Each one of these will produce a log file (`<prefix>_<module>.log`) and a subdirectory (`<prefix>_<module>`).

For this example, we will only look at the files within `GCF_000008865.2_summary`. The `<prefix>_summary` directory contains files that summarize the results that were generated across the pipeline. If you just want a list of the plasmids and viruses identified in your input, this is what you are looking for.

```
genomad_output
╰── GCF_000008865.2_summary
    ├── GCF_000008865.2_plasmid.fna
    ├── GCF_000008865.2_plasmid_genes.tsv
    ├── GCF_000008865.2_plasmid_proteins.faa
    ├── GCF_000008865.2_plasmid_summary.tsv
    ├── GCF_000008865.2_virus.fna
    ├── GCF_000008865.2_virus_genes.tsv
    ├── GCF_000008865.2_virus_proteins.faa
    ╰── GCF_000008865.2_virus_summary.tsv
```

First, let's look at `GCF_000008865.2_virus_summary.tsv`:

```
seq_name                               length   n_genes   genetic_code   virus_score   fdr   topology   coordinates       taxonomy
------------------------------------   ------   -------   ------------   -----------   ---   --------   ---------------   ------------------------------------------------------------
NC_002695.2|provirus_1245607_1309390   63784    90        11             0.9735        NA    Provirus   1245607-1309390   root;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes
NC_002695.2|provirus_891197_928364     37168    42        11             0.9704        NA    Provirus   891197-928364     root;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes
NC_002695.2|provirus_5041220_5079596   38377    53        11             0.9698        NA    Provirus   5041220-5079596   root;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes
…
```

This tabular file lists all the viruses that geNomad found in your input and gives you some convenient information about them. Here's what each column contains:

- `seq_name`: The identifier of the sequence in the input FASTA file. Proviruses will have the following name scheme: `<sequence_identifier>|provirus_<start_coordinate>_<end_coordinate>`.
- `length`: Length of the sequence (or the provirus, in the case of integrated viruses).
- `n_genes`: Number of genes encoded in the sequence.
- `genetic_code`: Predicted genetic code. Possible values are: 11 (standard code for Bacteria and Archaea), 4 (recoded TGA stop codon), or 15 (recoded TAG stop codon).
- `virus_score`: A measure of how confident geNomad is that the sequence is a virus. Sequences that have scores close to 1.0 are more likely to be viruses than the ones that have lower scores.
- `fdr`: The estimated false discovery rate (FDR) of the classification (that is, the expected proportion of false positives among the sequences up to this row). To estimate FDRs geNomad requires [score calibration](https://portal.nersc.gov/genomad/score_calibration.html), which is turned off by default. Therefore, this column will only contain `NA` values in this example.
- `topology`: Topology of the viral sequence. Possible values are: Linear, DTR (direct terminal repeats), ITR (inverted terminal repeats), or Provirus (viruses integrated in host genomes).
- `coordinates`: 1-indexed coordinates of the provirus region within host sequences. Will be `NA` for viruses that were not predicted to be integrated.
- `taxonomy`: Taxonomic assignment of the virus genome. Lineages follow the taxonomy contained in [ICTV's VMR number 19](https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/13426).

In our example, geNomad identified several proviruses integrated into the *E. coli* genome. They were all predicted to use the genetic code 11 and were assigned to the *Caudoviricetes* class, which contains all the tailed bacteriphages. Since they all have high scores, we can be confident that these are indeed viruses.

Another important file is `GCF_000008865.2_virus_genes.tsv`. During its execution, geNomad annotates the genes encoded by the input sequences using a database of chromosome, plasmid, and virus-specific markers. The `<prefix>_virus_genes.tsv` file summarizes the annotation of the genes encoded by the identified viruses.

```
gene                                     start    end      length   strand   gc_content   genetic_code   rbs_motif        marker              evalue       bitscore   uscg   taxid   taxname          annotation_accessions                 annotation_description
--------------------------------------   ------   ------   ------   ------   ----------   ------------   --------------   -----------------   ----------   --------   ----   -----   --------------   -----------------------------------   -------------------------
NC_002695.2|provirus_300073_325822_264   300073   301047   975      -1       0.476        11             AGGAG/GGAGG      NA                  NA           NA         0      1       NA               NA                                    NA
NC_002695.2|provirus_300073_325822_265   301423   301812   390      -1       0.392        11             GGA/GAG/AGG      NA                  NA           NA         0      1       NA               NA                                    NA
NC_002695.2|provirus_300073_325822_266   301940   302653   714      -1       0.461        11             None             GENOMAD.222303.VP   4.577e-09    59         0      2561    Caudoviricetes   PF06223.15;K10762                     Minor tail protein T
NC_002695.2|provirus_300073_325822_267   302754   302954   201      1        0.428        11             AGGAGG           GENOMAD.061471.VV   1.653e-15    71         0      2561    Caudoviricetes   PF09048.13;TIGR03339;K22302;COG1609   Cro
NC_002695.2|provirus_300073_325822_268   303073   303366   294      1        0.503        11             AGGA/GGAG/GAGG   GENOMAD.129061.VV   5.21e-33     123        0      2561    Caudoviricetes   PF05269.14;TIGR00721                  Bacteriophage CII protein
…
```

The columns in this file are:

- `gene`: Identifier of the gene (`<sequence_name>_<gene_number>`). Gene numbers start with 1 (first gene in the sequence). Because the genes in this example are encoded by prophages integrated in the middle of the host chromosome, they won't necessarily start with 1.
- `start`: 1-indexed start coordinate of the gene.
- `end`: 1-indexed end coordinate of the gene.
- `length`: Length of the gene locus (in base pairs).
- `strand`: Strand that encodes the gene. Can be 1 (direct strand) or -1 (reverse strand).
- `gc_content`: GC content of the gene locus.
- `genetic_code`: Predicted genetic code (see details in the explanation of the summary file).
- `rbs_motif`: Detected motif of the ribosome-binding site.
- `marker`: Best matching geNomad marker. If this gene doesn't match any markers, the value will be `NA`.
- `evalue`: E-value of the alignment between the protein encoded by the gene and the best matching geNomad marker.
- `bitscore`: Bitscore of the alignment between the protein encoded by the gene and the best matching geNomad marker.
- `uscg`: Whether the marker assigned to this gene corresponds to a universal single copy gene (UCSG). These genes are expected to be found in chromosomes and are rare in plasmids and viruses. Can be 1 (gene is USCG) or 0 (gene is not USCG).
- `taxid`: Taxonomic identifier of the marker assigned to this gene (you can ignore this as it is meant to be used internally by geNomad).
- `taxname`: Name of the taxon associated with the assigned geNomad marker. In this example, we can see that the annotated proteins are all characteristic of *Caudoviricetes* (which is why the provirus was assigned to this class).
- `annotation_accessions`: Some of the geNomad markers are functionally annotated. This column tells you which entries in Pfam, TIGRFAM, COG, and KEGG were assigned to the marker.
- `annotation_description`: A text describing the function assigned to the marker.

In the snippet above we can see the information of the first five genes encoded by `NC_002695.2|provirus_300073_325822`. The first two entried didn't match any geNomad marker. The following three were all assigned to protein families that are typical of tailed bacteriphages (such as the minor tail protein), which reassures us that these are indeed *Caudoviricetes*.

One important detail here is that the primary purpose of geNomad's markers is classification. They were designed to be specific to chromosomes, plasmids, or viruses, enabling the distinction of sequences belonging to these classes. Therefore, you should not expect that every single viral gene will be annotated with a geNomad marker. If you want to annotate the genes within your sequences as throughly as possible, you should use databases such as [Pfam](http://pfam.xfam.org/) or [EggNOG](http://eggnog5.embl.de/#/app/home).

The other two virus-related files within the summary directory are `GCF_000008865.2_virus.fna` and `CF_000008865.2_virus_proteins.faa`. These are FASTA files of the identified virus sequences and their proteins, respectively. Proviruses are automatically excised from the host sequence.

Enough with viruses. What about the plasmids?

As you would expect, the data pertaining to the identification of plasmids can be found in the `<prefix>_plasmid_summary.tsv`, `<prefix>_genes.tsv`, `<prefix>_plasmid.fna`, and `<prefix>_plasmid_proteins.faa` files. These are mostly very similar to their virus counterparts. The only difference is that `<prefix>_plasmid_summary.tsv` (shown below) doesn't have the virus-specific columns that are in `<prefix>_virus_summary.tsv` (`coordinates` and `topology`).

```
seq_name      length   n_genes   genetic_code   plasmid_score   fdr   topology
-----------   ------   -------   ------------   -------------   ---   --------
NC_002128.1   92721    88        11             0.9942          NA    Linear
NC_002127.1   3306     3         11             0.9913          NA    Linear
```
