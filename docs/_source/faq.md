# Frequently asked questions

Here you will find answers to some common questions about geNomad. If you can't find an answer to your problem here, please open an issue in the [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/apcamargo/genomad/).

## How can I speed up geNomad?

If you want to speed up the execution of geNomad, there are two options available:

- Disable the [neural network-based classification](nn_classification.md) using the `--disable-nn-classification` option, which will also disable [score aggregation](score_aggregation.md) and force geNomad to solely rely on the marker-based classifier.
- Decrease the sensitivity of the MMseqs2 search that assigns markers to genes with the `--sensitivity` parameter. This will make the [`annotate`](annotate-module) module faster, but will also decrease the number of genes assigned to markers.

Please note that both options may negatively impact geNomad’s classification performance.

## How can I get the taxonomy of all my sequences, regardless of their classification?

If you are only interested in performing taxonomic assignment of viral genomes and not classification of the sequences, you can run the [`annotate`](annotate-module) module alone using the following command:

```bash
genomad annotate input.fna genomad_output genomad_db
```

The taxonomic assignment of your sequences will be in the `genomad_output/input_annotate/input_taxonomy.tsv` file. Please note that the path may vary depending on the name of your input file, but the structure will remain consistent.

## Why is the execution is stuck at the *"Classifying sequences"* step of the `nn-classification` module?

During the [neural network-based classification](nn_classification.md) of sequences, the processor performs numerous computationally expensive matrix multiplications. Although modern processors can handle these computations quickly, some hardware, particularly older or laptop processors, can be quite slow during this step.

If you ran geNomad using the `end-to-end` command and your execution is stuck at the *"Classifying sequences"* step of the [`nn-classification`](nn-classification-module) module, you can disable the neural network classifier by using the `--disable-nn-classification` option. By doing this, geNomad will rely solely on marker-based classification and disable [score aggregation](score_aggregation.md).

## Why execution is failing at the *"Annotating proteins with MMseqs2…"* step of the `annotate` module?

If geNomad fails during the *"Annotating proteins with MMseqs2 and geNomad database"* step, it is likely that your computer is running out of memory. This can be resolved by using the[`--splits`](notes-about-parameters) parameter, which splits the marker database into smaller chuncks and searches each of them independently.
