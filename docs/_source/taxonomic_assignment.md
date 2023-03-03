# Taxonomic assignment of virus genomes

geNomad is able to assign viral genomes to taxonomic lineages defined in [ICTV's VMR number 19](https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/13426). This is allowed by the 85,315 markers that are associated with viral taxa that cover most of the lineages recognized by the [ICTV](https://talk.ictvonline.org/).

## Taxonomic assignment algorithm

To assign viral sequences to specific taxa, geNomad first aligns the genes encoded by the sequences to a set of 227,897 markers. The genes that produce significant matches are then assigned to a marker, which might contain taxonomic information (colored genes below).

```{image} _static/figures/taxonomic_assignment_1.svg
:width: 420
:class: no-scaled-link
:align: center
```

Each gene is subsequently classified based on the taxonomic lineage of the assigned marker. Different genes within the sequence might be assigned to different lineages.

```{image} _static/figures/taxonomic_assignment_2.svg
:width: 1000
:align: center
```

To establish a single sequence-level taxonomy, weights are computed for each taxon included among the gene-level assignments, as well as their parent taxa (up to the root of the taxonomy) by summing the bitscores obtained from the alignments with marker profiles.

```{image} _static/figures/taxonomic_assignment_3.svg
:width: 1000
:align: center
```

The taxonomy of the sequence is determined as the most specific taxon that is supported by at least 50% of the total weight (sum of the bitscores of all genes with taxonomy). In the example above, both families (*Autographiviridae* and *Zobellviridae*) failed to reach 50% consensus, so the genome is assigned to the class *Caudoviricetes*.

**Final classification:** *Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes*.