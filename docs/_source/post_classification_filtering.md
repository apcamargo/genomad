# Post-classification filtering

To identify plasmids and viruses in a set of sequences, geNomad uses classification models that classify sequences into chromosomes, plasmids, or viruses. However, relying solely on the model's outputs can be problematic. For example, short sequences (less than 2,500 bp, for example) are difficult to classify and are often misclassified, requiring extra caution. Additionally, different users may have different expectations regarding results; some may want as many predictions as possible, while others may prefer only the most strongly supported predictions.

To address these challenges, geNomad applies a series of filters to the classification results. This ensures that users receive reliable predictions while allowing them to control the level of stringency of the filtering process.

## Filters

The following filters are available to users to generate the final lists of plasmids and viruses:

- **`min-score`:** Minimum score to flag a sequence as virus or plasmid.
- **`max-fdr`:** Maximum accepted false discovery rate. This option will be ignored if the scores were not calibrated.
- **`min-number-genes`:** The minimum number of genes a sequence must encode to be considered for classification as a plasmid or virus.
- **`min-plasmid-marker-enrichment`:** Minimum allowed value for the plasmid marker enrichment score.
- **`min-virus-marker-enrichment`:** Minimum allowed value for the virus marker enrichment score.
- **`min-plasmid-hallmarks`:** Minimum number of plasmid hallmarks in the identified plasmids.
- **`min-plasmid-hallmarks-short-seqs`:** Minimum number of plasmid hallmarks in the identified plasmids that are shorter than 2,500 bp.
- **`min-virus-hallmarks`:** Minimum number of virus hallmarks in the identified viruses.
- **`min-virus-hallmarks-short-seqs`:** Minimum number of virus hallmarks in the identified viruses that are shorter than 2,500 bp.
- **`max-uscg`:** Maximum allowed number of universal single copy genes (USCGs) in a virus or a plasmid.

```{admonition} The marker enrichment filters
:class: tip
The `min-plasmid-marker-enrichment` and `min-virus-marker-enrichment` filters rely on the computation of the marker enrichment of each sequence, which represents the total enrichment of plasmid or virus markers in it. Essentially, sequences with multiple markers will have higher values than the ones that encode few or no markers.

$$\textrm{Plasmid marker enrichment} = \sum _{i=1}^ne^{\textrm{P SPM}_i}-e^{\textrm{C SPM}_i+\textrm{V SPM}_i}\:$$

$$\textrm{Virus marker enrichment} = \sum _{i=1}^ne^{\textrm{V SPM}_i}-e^{\textrm{C SPM}_i+\textrm{P SPM}_i}\:$$
```

(default-parameters-and-presets)=
## Default parameters and presets

Given the large number of available filters, geNomad also allows users to use presets that can disable all filters (`--relaxed`) or make the filtering process more aggressive (`--conservative`).

The values used to filter predictions when executing geNomad with default parameters or one of the presets are the following:

| Filter                             | Default | Relaxed | Conservative |
|:-----------------------------------|--------:|--------:|-------------:|
| `min-score`                        |    0.70 |    0.00 |         0.80 |
| `max-fdr`                          |    0.10 |    1.00 |         0.05 |
| `min-number-genes`                 |       1 |       0 |            1 |
| `min-plasmid-marker-enrichment`    |    0.10 | -100.00 |         1.50 |
| `min-virus-marker-enrichment`      |    0.00 | -100.00 |         1.50 |
| `min-plasmid-hallmarks`            |       0 |       0 |            1 |
| `min-plasmid-hallmarks-short-seqs` |       1 |       0 |            1 |
| `min-virus-hallmarks`              |       0 |       0 |            1 |
| `min-virus-hallmarks-short-seqs`   |       1 |       0 |            1 |
| `max-uscg`                         |       4 |     100 |            2 |

```{admonition} Post-classification filtering in the absence of gene annotation
:class: tip
When the `annotate` module is not executed, most filters are disabled, since they rely on gene annotation. In such cases, only `min-score` and `max-fdr` are used to filter the classification results.
```
