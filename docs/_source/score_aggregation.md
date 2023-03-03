# Hybrid classification framework and score aggregation

geNomad employs a hybrid approach to plasmid and virus identification that combines a neural [network-based classifier](nn_classification.md), that dispenses alignments to reference databases, and a marker-based classifier, that classifies sequences based on the presence of [protein markers](marker_features.md) that are informative for classification.

```{image} _static/figures/score_aggregation_1.svg
:width: 670
:class: no-scaled-link
:align: center
```

To improve classification performance, geNomad capitalizes on the strengths of both models by aggregating their outputs into a single classification.

## Attention-based score aggregation

The neural network-based and the marker-based classifiers use distinct and often complementary approaches to classify input sequences. By aggregating their outputs, geNomad can take advantage of both approaches and provide a better classification performance. This is achieved through an attention mechanism that consists of a linear model that weighs the branches based on the frequency of chromosome, plasmid, and virus markers in the input sequence.

```{image} _static/figures/score_aggregation_2.svg
:width: 380
:class: no-scaled-link
:align: center
```

The attention mechanism works in such a way that the contribution of the marker branch goes higher as the fraction of genes that are assigned to markers increases. Essentially, the branch aggregation gives more weight to the marker branch when it is more informative (i.e., when most of the genes encoded by the input sequence are assigned to markers) and relies more on the sequence branch when marker information is scarce.

```{image} _static/figures/score_aggregation_3.svg
:width: 320
:class: no-scaled-link
:align: center
```
