# Provirus identification

Temperate phages can integrate into host genomes and form proviruses, which can greatly affect the host metabolism and ecology. geNomad is able to estimate the coordinates of integrated proviruses within host sequences, providing valuable information to the user.

## Provirus demarcation algorithm

geNomad's marker dataset provides information on how specific each marker is to viruses or hosts (which you can read more about [here](marker_features.md)). To identify putative proviruses, geNomad identifies sequences that contain regions of virus-specific markers that surrounded by host-specific genes.

Therefore, the first step in the process is assigning the genes encoded by a given sequence to geNomad markers,

```{image} _static/figures/provirus_1.svg
:width: 480
:class: no-scaled-link
:align: center
```

To demarcate regions that possibly correspond to proviruses, geNomad employs a conditional random field (CRF) model. This model takes into account the chromosome and virus SPM values of genes annotated with geNomad markers as input and uses contextual information to calculate the conditional probability of a sequence of states (either chromosome or provirus) in the sequence.

```{image} _static/figures/provirus_2.svg
:width: 480
:class: no-scaled-link
:align: center
```

As a result, the CRF model provides, for each gene, a a probability of it belonging to a provirus.

```{image} _static/figures/provirus_3.svg
:width: 480
:class: no-scaled-link
:align: center
```

Genes are then assigned to their most likely states, forming provirus islands that represent regions that are enriched in virus markers.

```{image} _static/figures/provirus_4.svg
:width: 480
:class: no-scaled-link
:align: center
```

To prevent having proviruses split into multiple islands due to incomplete marker coverage, provirus islands that are separated by short gene arrays (less than 6 genes or 2 chromosome markers) are merged. Next, the total viral enrichment of each island is computed as follows:

$$\textrm{Total viral enrichment} = \sum _{i=1}^ne^{\textrm{V SPM}_i}-e^{\textrm{C SPM}_i}\:$$

These values represent the total viral enrichment of the islands, taking into account all of the genes within them. Islands with several virus-specific markers will have higher marker enrichment, while islands with few virus-specific genes will have low marker enrichment.

```{image} _static/figures/provirus_5.svg
:width: 480
:class: no-scaled-link
:align: center
```

After identifying viral islands, geNomad refines their boundaries using tRNA and integrase annotation. This is because tRNAs and integrases are often located next to integrated elements due to site-specific recombination dynamics. Therefore, geNomad extends the provirus boundaries to the neighboring tRNAs (identified using [`ARAGORN`](http://www.ansikte.se/ARAGORN/)) and/or integrases (identified with [`MMseqs2`](https://github.com/soedinglab/MMseqs2/), using a set of 16 profiles of site-specific tyrosine integrases). Finally, islands with low marker enrichment are filtered out as they are usually not *bona fide* proviruses.

```{image} _static/figures/provirus_6.svg
:width: 480
:class: no-scaled-link
:align: center
```
