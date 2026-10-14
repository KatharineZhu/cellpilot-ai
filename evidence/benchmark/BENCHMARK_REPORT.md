# Independent annotation benchmark

**Status:** completed on real public FACS-enriched 10x datasets

## Measured strict metrics

| Metric | Value |
|---|---:|
| `n_cells` | 2500 |
| `coverage` | 0.965600 |
| `strict_accuracy` | 0.932800 |
| `strict_macro_f1` | 0.947814 |
| `strict_balanced_accuracy` | 0.932800 |
| `strict_mcc` | 0.918448 |
| `runtime_seconds` | 1.057777 |

## Reference-label design

Each reference class is the named population in an official 10x FACS-enriched dataset. Labels are held in an external Series and are absent from the AnnData object supplied to CellTypist.

## Limitations

- FACS enrichment labels are experimental population labels, not perfect single-cell ground truth; reported purity ranges from 92% to 100%.
- Each broad class is a separate 10x library (and may be a separate donor), so biological class is confounded with library/donor and technical effects.
- Inputs are Cell Ranger filtered gene-barcode matrices after cell calling, not FASTQ files or unfiltered droplets.
- The benchmark covers only five broad blood lineages and does not establish fine-subtype performance or out-of-reference rejection performance.
- The frozen CellTypist model cites Stephenson et al. 2021 rather than these 10x panels; this reduces known study-level leakage but does not prove the absence of every indirect training-data relationship.
- Cell-level metrics do not estimate donor-to-donor generalization because there are no independent biological replicates within each class panel.

Unmapped predictions count as errors in every strict metric. These results must not be described as donor-level generalization or fine-subtype accuracy.
