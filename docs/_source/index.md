---
hide-toc: true
---

# geNomad

geNomad is a tool that identifies virus and plasmid genomes from nucleotide sequences. It provides state-of-the-art classification performance and can be used to quickly find mobile genetic elements from genomes, metagenomes, or metatranscriptomes.

::::{grid}
:gutter: 2

:::{grid-item-card}
:columns: 12 12 4 4
**Speed**
^^^
geNomad is significantly faster than similar tools and can be used to process large datasets.
:::

:::{grid-item-card}
:columns: 12 12 4 4
**Taxonomic assignment**
^^^
The identified viruses are assigned to taxonomic lineages that follow the latest [ICTV](https://talk.ictvonline.org/) taxonomy release.
:::

:::{grid-item-card}
:columns: 12 12 4 4
**Functional annotation**
^^^
Genes encoded by viruses and plasmids are functionally annotated using geNomad's marker database.
:::

::::

## {octicon}`rocket;0.85em` Get started

To start using geNomad, read the installation and quickstart guides below. In case you want to learn more details about how geNomad works, visit the pages listed in the sidebar.

:::{card} Installation
:link: installation
:link-type: doc
Instructions on how to install geNomad in your computer or server.
:::

:::{card} Quickstart
:link: quickstart
:link-type: doc
Learn how to run geNomad and interpret its results.
:::

## {octicon}`browser;0.85em` Web app

```{image} _static/figures/nmdc_edge_logo.png
:width: 210
:class: no-scaled-link
:align: center
:target: https://nmdc-edge.org/virus_plasmid/workflow
```

geNomad is also available as a web app in the [NMDC EDGE](https://nmdc-edge.org/virus_plasmid/workflow) platform. There you can upload your sequence data, visualize the results in your browser, and download the data to your computer.

## {octicon}`bookmark;0.85em` Citing geNomad

If you use geNomad in your work, please consider citing its manuscript:

:::{card}
:link: https://www.nature.com/articles/s41587-023-01953-y

**Identification of mobile genetic elements with geNomad**
Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C. — *Nature Biotechnology* (2023), DOI: 10.1038/s41587-023-01953-y.
:::

## {octicon}`question;0.85em` Ask a question or report a bug

If you want to ask a question about geNomad or report a problem you had with it, please create an issue in the [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/apcamargo/genomad/).

```{toctree}
:hidden:

self
```

```{toctree}
:caption: Using geNomad
:hidden:

installation
quickstart
pipeline
faq
```

```{toctree}
:caption: Theory
:hidden:

score_aggregation
marker_features
nn_classification
provirus_identification
taxonomic_assignment
score_calibration
post_classification_filtering
```
