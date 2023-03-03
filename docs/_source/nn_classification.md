# Neural network-based classification

Neural networks can be used to perform classification directly from nucleotide sequences, dispensing explicit alignments to reference genomes or proteins. This is possible because some neural network architectures are capable of learning discriminative sequence motifs that are informative for classification. In geNomad, a neural network-based classifier is used in conjunction with a marker-based classifier to improve classification of sequences with few or no genes assigned to markers.

## geNomad's neural network classifier

To identify sequences of plasmids and viruses without explicit alignment to a reference database, geNomad processes inputs using a neural network model that is able to classify the sequences from their nucleotide makeup alone. To accomplish this, the input sequences are first converted to a numerical format by tokenizing them into arrays of 4-mer words, which are then one-hot-encoded, creating binary 256-dimensional matrices that reflect the presence of specific 4-mers (rows) across different positions within the sequence (columns).

```{image} _static/figures/nn_classification_1.svg
:width: 680
:class: no-scaled-link
:align: center
```

These matrices are then passed to an encoder, which generates vector representations of the sequences in a high-dimensional embedding space. In this space, sequence representations from the same class (chromosome, plasmid, or virus) will be more similar compared to sequence representations from different classes. The resulting representations are subsequently fed to a dense neural network that produces the classification scores.

```{image} _static/figures/nn_classification_2.svg
:width: 580
:class: no-scaled-link
:align: center
```
## The IGLOO architecture

To generate vector representations of the inputs, geNomad employs an encoder based on the [IGLOO](https://arxiv.org/abs/1807.03402) architecture, which is able to extract patterns that are useful for classification from the sequence data and encode them into an embedding space. The IGLOO encoder begins processing one-hot-encoded matrices by applying convolutional filters to generate sequence feature maps. To capture relationships between non-contiguous parts of the sequence, IGLOO generates patches that contain slices taken from random positions within the sequence.

```{image} _static/figures/nn_classification_3.svg
:width: 450
:class: no-scaled-link
:align: center
```

These patches are integrated into a self-attention mechanism that weights different parts of the feature map and uses the long-range dependencies encoded in the patches to derive the final sequence representation.

```{image} _static/figures/nn_classification_4.svg
:width: 630
:class: no-scaled-link
:align: center
```
