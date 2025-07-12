# geNomad

geNomad: Identification of mobile genetic elements

## Features

geNomad's primary goal is to identify viruses and plasmids in sequencing data (isolates, metagenomes, and metatranscriptomes). It also provides a couple of additional features that can help you in your analysis:

- Taxonomic assignment of viral genomes.
- Identification of viruses integrated in host genomes (proviruses).
- Functional annotation of proteins.

## Documentation

For installation instructions, information about how geNomad works, and a detailed explanation of how to execute it, please check the full documentation: https://portal.nersc.gov/genomad/

## Web app

geNomad is available as a web application on the [Galaxy](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/ufz/genomad_end_to_end/genomad_end_to_end/1.11.1+galaxy0) and [NMDC EDGE](https://nmdc-edge.org/virus_plasmid/workflow) platforms. There, you can upload your sequence data, view the results directly in your browser, and download the outputs for further analysis.


## Citing geNomad

If you use geNomad in your work, please consider citing its manuscript:

> [**Identification of mobile genetic elements with geNomad**](https://www.nature.com/articles/s41587-023-01953-y)
>
> Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C. — *Nature Biotechnology* (2023), DOI: 10.1038/s41587-023-01953-y
        
        
        
        .

## Quick start

We recommend users to read the [documentation](https://portal.nersc.gov/genomad/) before starting to use geNomad. If you are in a rush, however, you can follow this quick step-by-step example.

### Installation

First, you need to install geNomad. There's a couple of ways to do that (read more about it in the [documentation](https://portal.nersc.gov/genomad/installation.html)), but two convinient options are using [Pixi](https://pixi.sh/) or [Mamba](https://mamba.readthedocs.io/). Both of them will handle the installation of all dependencies for you.

Pixi allows you to install geNomad as a globally available command for easy execution.

```
pixi global install -c conda-forge -c bioconda genomad
```

With Mamba, you will create an environment for geNomad and activate it before being able to use it.

```
# Create an environment for geNomad
mamba create -n genomad -c conda-forge -c bioconda genomad
# Activate the geNomad environment
mamba activate genomad
```

Another option is to use geNomad through Docker.

```
# Pull the image
docker pull antoniopcamargo/genomad
# Run the image
docker run --rm -ti -v "$(pwd):/app" antoniopcamargo/genomad
```

### Downloading the database

geNomad depends on a database that contains the profiles of the markers that are used to classify sequences, their taxonomic information, their functional annotation, etc. So, you should first download the database to your current directory:

```
genomad download-database .
```

The database will be contained within the `genomad_db` directory.

If you prefer, you can also download the database from [Zenodo](https://zenodo.org/records/14886553) and extract it manually.

### Executing geNomad

Now you are ready to go! geNomad works by executing a series of modules sequentially (you can find more information about this in the [pipeline documentation](https://portal.nersc.gov/genomad/pipeline.html)), but we provide a convenient `end-to-end` command that will execute the entire pipeline for you in one go.

In this example, we will use an *Klebsiella pneumoniae* genome ([GCF_009025895.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009025895.1/)) as input. You can use any FASTA file containing nucleotide sequences as input. geNomad will work for isolate genomes, metagenomes, and metatranscriptomes.

The command to execute geNomad is structured like this:

```
genomad end-to-end [OPTIONS] INPUT OUTPUT DATABASE
```

So, to run the full geNomad pipeline (`end-to-end` command), taking a nucleotide FASTA file (`GCF_009025895.1.fna.gz`) and the database (`genomad_db`) as input, we will execute the following command:

```
genomad end-to-end --cleanup --splits 8 GCF_009025895.1.fna.gz genomad_output genomad_db
```

The results will be written inside the `genomad_output` directory.

Three important details about the command above:

- The `--cleanup` option was used to force geNomad to delete intermediate files that were generated during the execution. This will save you some storage space.
- The `--splits 8` parameter was used here to make it possible to run this example in a notebook. geNomad searches a big database of protein profiles that take up a lot of space in memory. To prevent the execution from failing due to insufficient memory, we can use the `--splits` parameter to split the search into chuncks. If you are running geNomad in a big server you might not need to split your search, increasing the execution speed.
- Note that the input FASTA file that I used as input was compressed. This is possible because geNomad supports input files compressed as `.gz`, `.bz2`, or `.xz`.

> [!NOTE]
> By default, geNomad applies a series of post-classification filters to remove likely false positives. For example, sequences are required to have a plasmid or virus score of at least 0.7 and sequences shorter than 2,500 bp are required to encode at least one hallmark gene. If you want to disable the post-classification filters, add the `--relaxed` flag to your command. On the other hand, if you want to be very conservative with your classification, you may use the `--conservative` flag. This will make the post-classification filters more aggressive, preventing sequences without strong support from being classified as plasmid or virus. You can check out the default, relaxed, and conservative post-classification filters [here](https://portal.nersc.gov/genomad/post_classification_filtering.html#default-parameters-and-presets).

### Understanding the outputs

In this example, the results of geNomad's analysis will be written to the `genomad_output` directory, which will look like this:

```
genomad_output
├── GCF_009025895.1_aggregated_classification
├── GCF_009025895.1_aggregated_classification.log
├── GCF_009025895.1_annotate
├── GCF_009025895.1_annotate.log
├── GCF_009025895.1_find_proviruses
├── GCF_009025895.1_find_proviruses.log
├── GCF_009025895.1_marker_classification
├── GCF_009025895.1_marker_classification.log
├── GCF_009025895.1_nn_classification
├── GCF_009025895.1_nn_classification.log
├── GCF_009025895.1_summary
╰── GCF_009025895.1_summary.log
```

As mentioned above, geNomad works by executing several modules sequentially. Each one of these will produce a log file (`<prefix>_<module>.log`) and a subdirectory (`<prefix>_<module>`).

For this example, we will only look at the files within `GCF_009025895.1_summary`. The `<prefix>_summary` directory contains files that summarize the results that were generated across the pipeline. If you just want a list of the plasmids and viruses identified in your input, this is what you are looking for.

```
genomad_output
╰── GCF_009025895.1_summary
    ├── GCF_009025895.1_plasmid.fna
    ├── GCF_009025895.1_plasmid_genes.tsv
    ├── GCF_009025895.1_plasmid_proteins.faa
    ├── GCF_009025895.1_plasmid_summary.tsv
    ├── GCF_009025895.1_summary.json
    ├── GCF_009025895.1_virus.fna
    ├── GCF_009025895.1_virus_genes.tsv
    ├── GCF_009025895.1_virus_proteins.faa
    ╰── GCF_009025895.1_virus_summary.tsv
```

First, let's look at `GCF_009025895.1_virus_summary.tsv`:

```
seq_name                                 length   topology              coordinates       n_genes   genetic_code   virus_score   fdr   n_hallmarks   marker_enrichment   taxonomy
--------------------------------------   ------   -------------------   ---------------   -------   ------------   -----------   ---   -----------   -----------------   -----------------------------------------------------------------
NZ_CP045015.1|provirus_3855947_3906705   50759    Provirus              3855947-3906705   79        11             0.9772        NA    16            73.7974             Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;
NZ_CP045015.1|provirus_2885031_2934610   49580    Provirus              2885031-2934610   70        11             0.9769        NA    13            73.2757             Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;
NZ_CP045018.1                            51887    No terminal repeats   NA                57        11             0.9760        NA    14            65.4720             Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;
…
```

This tabular file lists all the viruses that geNomad found in your input and gives you some convenient information about them. Here's what each column contains:

- `seq_name`: The identifier of the sequence in the input FASTA file. Proviruses will have the following name scheme: `<sequence_identifier>|provirus_<start_coordinate>_<end_coordinate>`.
- `length`: Length of the sequence (or the provirus, in the case of integrated viruses).
- `topology`: Topology of the viral sequence. Possible values are: `No terminal repeats`, `DTR` (direct terminal repeats), `ITR` (inverted terminal repeats), or `Provirus` (viruses integrated in host genomes).
- `coordinates`: 1-indexed coordinates of the provirus region within host sequences. Will be `NA` for viruses that were not predicted to be integrated.
- `n_genes`: Number of genes encoded in the sequence.
- `genetic_code`: Predicted genetic code. Possible values are: 11 (standard code for Bacteria and Archaea), 4 (recoded TGA stop codon), or 15 (recoded TAG stop codon).
- `virus_score`: A measure of how confident geNomad is that the sequence is a virus. Sequences that have scores close to 1.0 are more likely to be viruses than the ones that have lower scores.
- `fdr`: The estimated false discovery rate (FDR) of the classification (that is, the expected proportion of false positives among the sequences up to this row). To estimate FDRs geNomad requires [score calibration](https://portal.nersc.gov/genomad/score_calibration.html), which is turned off by default. Therefore, this column will only contain `NA` values in this example.
- `n_hallmarks`: Number of genes that matched a hallmark geNomad marker. Hallmarks are genes that were previously associated with viral function and their presence is a strong indicative that the sequence is indeed a virus.
- `marker_enrichment`: A score that represents the total enrichment of viral markers in the sequence. The value goes as the number of virus markers in the sequence increases, so sequences with multiple markers will have higher score. Chromosome and plasmid markers will reduce the score.
- `taxonomy`: Taxonomic assignment of the virus genome. Lineages follow the taxonomy contained in [ICTV's Taxonomy Release MSL39](https://ictv.global/sites/default/files/VMR/VMR_MSL39.v4_20241106.xlsx). Viruses can be taxonomically assigned up to the family level, but not to specific genera or species within that family. The taxonomy is presented with a fixed number of fields (corresponding to taxonomic ranks) separated by semicolons, with empty fields left blank.

In our example, geNomad identified several proviruses integrated into the *K. pneumoniae* genome and one extrachromosomal phage. Since they all have high scores and marker enrichment, we can be confident that these are indeed viruses. They were all predicted to use the genetic code 11 and were assigned to the *Caudoviricetes* class, which contains all the tailed bacteriphages. In the `taxonomy` field for these viruses, after *Caudoviricetes*, there are two consecutive semicolons because geNomad could only assign them to the class level, leaving the order and family ranks empty.

Another important file is `GCF_009025895.1_virus_genes.tsv`. During its execution, geNomad annotates the genes encoded by the input sequences using a database of chromosome, plasmid, and virus-specific markers. The `<prefix>_virus_genes.tsv` file summarizes the annotation of the genes encoded by the identified viruses.

```
gene              start   end     length   strand   gc_content   genetic_code   rbs_motif     marker              evalue       bitscore   uscg   plasmid_hallmark   virus_hallmark   taxid   taxname          annotation_conjscan   annotation_amr   annotation_accessions              annotation_description
---------------   -----   -----   ------   ------   ----------   ------------   -----------   -----------------   ----------   --------   ----   ----------------   --------------   -----   --------------   -------------------   --------------   --------------------------------   --------------------------------------------------------------------------------------
NZ_CP045018.1_1   1       399     399      1        0.536        11             None          GENOMAD.108715.VP   2.675e-31    120        0      0                  1                1246    Caudoviricetes   NA                    NA               PF05100;COG4672;TIGR01600          Phage minor tail protein L
NZ_CP045018.1_2   401     1111    711      1        0.568        11             AGGAG         GENOMAD.168265.VP   1.523e-39    149        0      0                  0                1246    Caudoviricetes   NA                    NA               PF14464;COG1310;K21140;TIGR02256   Proteasome lid subunit RPN8/RPN11, contains Jab1/MPN domain metalloenzyme (JAMM) motif
NZ_CP045018.1_3   1143    1493    351      1        0.382        11             AGGAG         GENOMAD.147875.VV   7.818e-13    66         0      0                  0                1246    Caudoviricetes   NA                    NA               COG5633;TIGR03066                  NA
NZ_CP045018.1_4   1509    2120    612      1        0.477        11             GGA/GAG/AGG   GENOMAD.143103.VP   2.238e-48    173        0      0                  1                1246    Caudoviricetes   NA                    NA               PF06805;COG4723;TIGR01687          Phage-related protein, tail component
NZ_CP045018.1_5   2183    13516   11334    1        0.566        11             None          GENOMAD.159864.VP   7.104e-262   901        0      0                  0                1246    Caudoviricetes   NA                    NA               PF12421;PF09327                    Fibronectin type III protein
NZ_CP045018.1_6   13585   15084   1500     1        0.550        11             AGGAG         GENOMAD.195756.VP   1.695e-13    76         0      0                  0                1246    Caudoviricetes   NA                    NA               NA                                 NA
NZ_CP045018.1_7   15163   16128   966      -1       0.469        11             GGAGG         NA                  NA           NA         0      0                  0                1       NA               NA                    NA               NA                                 NA
…
```

The columns in this file are:

- `gene`: Identifier of the gene (`<sequence_name>_<gene_number>`). Usually, gene numbers start with 1 (first gene in the sequence). However, genes encoded by prophages integrated in the middle of the host chromosome may start with a different number, depending on it's position within the chromosome.
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
- `uscg`: Whether the marker assigned to this gene corresponds to a universal single-copy gene (UCSG, as defined in [BUSCO v5](https://busco.ezlab.org/)). These genes are expected to be found in chromosomes and are rare in plasmids and viruses. Can be 1 (gene is USCG) or 0 (gene is not USCG).
- `plasmid_hallmark`: Whether the marker assigned to this gene represents a plasmid hallmark.
- `virus_hallmark`: Whether the marker assigned to this gene represents a virus hallmark.
- `taxid`: Taxonomic identifier of the marker assigned to this gene (you can ignore this as it is meant to be used internally by geNomad).
- `taxname`: Name of the taxon associated with the assigned geNomad marker. In this example, we can see that the annotated proteins are all characteristic of *Caudoviricetes* (which is why the provirus was assigned to this class).
- `annotation_conjscan`: If the marker that matched the gene is a conjugation-related gene (as defined in [CONJscan](https://link.springer.com/protocol/10.1007/978-1-4939-9877-7_19)) this field will show which CONJscan acession was assigned to the marker.
- `annotation_amr`: If the marker that matched the gene was annotated with an antimicrobial resistance (AMR) function (as defined in [NCBIfam-AMRFinder](https://www.ncbi.nlm.nih.gov/pathogens/hmm/)), this field will show which NCBIfam acession was assigned to the marker.
- `annotation_accessions`: Some of the geNomad markers are functionally annotated. This column tells you which entries in Pfam, TIGRFAM, COG, and KEGG were assigned to the marker.
- `annotation_description`: A text describing the function assigned to the marker.

In the example above we can see the information of the first seven genes encoded by `NZ_CP045018.1`. The last entry didn't match any geNomad marker. The first six were all assigned to protein families, some of which are typical of tailed bacteriphages (such as the minor tail protein), reassuring us that these are indeed *Caudoviricetes*.

One important detail here is that the primary purpose of geNomad's markers is classification. They were designed to be specific to chromosomes, plasmids, or viruses, enabling the distinction of sequences belonging to these classes. Therefore, you should not expect that every single viral gene will be annotated with a geNomad marker. If you want to annotate the genes within your sequences as throughly as possible, you should use databases such as [Pfam](http://pfam.xfam.org/) or [COG](https://www.ncbi.nlm.nih.gov/research/cog/).

The other two virus-related files within the summary directory are `GCF_009025895.1_virus.fna` and `GCF_009025895.1_virus_proteins.faa`. These are FASTA files of the identified virus sequences and their proteins, respectively. Proviruses are automatically excised from the host sequence.

Moving on to plasmids, the data related to their identification can be found in the `<prefix>_plasmid_summary.tsv`, `<prefix>_genes.tsv`, `<prefix>_plasmid.fna`, and `<prefix>_plasmid_proteins.faa` files. These are mostly very similar to their virus counterparts. The differences in `<prefix>_plasmid_summary.tsv` (shown below) are the following:

- Virus-specific columns that are in `<prefix>_virus_summary.tsv` (`coordinates` and `taxonomy`) are not present.
- The `conjugation_genes` column lists genes that might be involved in conjugation. It's important to note that the presence of such genes is not sufficient to tell whether a given plasmid is conjugative or mobilizible. If you are interested in identifying conjugative plasmids, we recommend you to analyze the plasmids you identified using geNomad with [CONJscan](https://link.springer.com/protocol/10.1007/978-1-4939-9877-7_19).
- The `amr_genes` column lists genes annotated with antimicrobial resistance function. You can check the specific functions associated with each accession in [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/hmm/) website.

```
seq_name        length   topology              n_genes   genetic_code   plasmid_score   fdr   n_hallmarks   marker_enrichment   conjugation_genes                                                                                       amr_genes
-------------   ------   -------------------   -------   ------------   -------------   ---   -----------   -----------------   -----------------------------------------------------------------------------------------------------   -----------------
NZ_CP045020.1   28729    No terminal repeats   36        11             0.9954          NA    6             26.4290             F_traE                                                                                                  NA
NZ_CP045022.1   50635    No terminal repeats   61        11             0.9946          NA    9             44.8458             T_virB1;T_virB3;virb4;T_virB5;T_virB6;T_virB8;T_virB9                                                   NA
NZ_CP045019.1   44850    No terminal repeats   52        11             0.9943          NA    4             27.8509             F_traE                                                                                                  NA
NZ_CP045016.1   82240    No terminal repeats   110       11             0.9935          NA    11            36.6048             T_virB8;T_virB9;F_traF;F_traH;F_traG;T_virB1                                                            NF000270;NF012171
NZ_CP045021.1   5251     No terminal repeats   7         11             0.9932          NA    2             3.1408              MOBP1                                                                                                   NA
NZ_CP045017.1   61331    No terminal repeats   76        11             0.9929          NA    16            35.2570             I_trbB;I_trbA;MOBP1;I_traI;I_traK;I_traL;I_traN;I_traO;I_traP;I_traQ;I_traR;traU;I_traW;I_traY;F_traE   NA
```
