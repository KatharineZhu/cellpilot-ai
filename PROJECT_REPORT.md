# CellPilot AI Project Report

Report date: 2026-10-14
Evidence snapshot: `results/pbmc3k/verified_run`

## Executive summary

CellPilot AI completes a real PBMC 3k analysis from a checksum-verified 10x
Genomics count archive through QC, normalization, PCA, UMAP, Leiden clustering,
marker ranking, CellTypist annotation, deterministic evidence-bounded
interpretation, and a five-class external benchmark. The saved run contains real
tables, plots, a processed AnnData object, reference-data manifests, stage
records, and an artifact-integrity manifest.

The evidence supports calling this a working alpha research workflow. It does
not support calling it production-ready, clinically validated, or biologically
definitive. A final independent scientific audit has not been performed.

## Verification basis

This report uses only on-disk evidence:

- PBMC and five benchmark source manifests under `data/`;
- the frozen CellTypist model metadata;
- analysis, annotation, interpretation, and benchmark artifacts beneath
  `results/pbmc3k/verified_run`;
- a fresh local `pytest -q` execution on 2026-10-14;
- a fresh local manifest verification on 2026-10-14.

No metric is inferred from a test fixture, cluster-consensus table, or planned
protocol.

## Environment

The verified manifest records CPython 3.12.5 on arm64 macOS with Scanpy 1.12.4,
AnnData 0.13.4, CellTypist 1.7.1, NumPy 2.5.4, pandas 3.0.6, SciPy 1.18.1,
scikit-learn 1.7.2, Streamlit 1.65.0, Pydantic 2.14.1, and OpenAI 2.54.0.

Python 3.14.4 was available on the host but was not used. Compatibility of this
stack with Python 3.14 remains untested.

## PBMC 3k data provenance

The main input is the public 10x Genomics PBMC 3k Cell Ranger-filtered
gene-barcode UMI matrix, not raw droplets or FASTQ files.

| Property | Recorded value |
| --- | --- |
| Source | Official 10x Genomics archive |
| License recorded from source page | CC BY 4.0 |
| Project-frozen and observed SHA-256 | `847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5` |
| Archive size | 7,621,991 bytes |
| Matrix shape | 2,700 cells × 32,738 genes |
| Nonzero count entries | 2,286,884 |
| Maximum UMI count | 419 |

The loader verifies finite, non-negative, integer-like counts before analysis.

## Analysis methods

### Quality control

The prespecified configuration retains cells with 200–2,500 detected genes and
less than 5% mitochondrial counts, then genes detected in at least three cells.
No total-count threshold is applied. These demonstration-specific thresholds
were not optimized against benchmark outcomes and are not universal guidance.

### Normalization and feature selection

The workflow copies filtered raw counts to `layers["counts"]`, normalizes each
cell to 10,000, and applies `log1p` once. It stores the full-gene log-normalized
matrix in `.raw`. Scaling occurs only in a temporary HVG subset for PCA, leaving
the main expression matrix non-negative for CellTypist. The run selected 2,000
HVGs with Scanpy's `seurat` flavor.

### Embedding, clustering, and markers

PCA used 40 components; the neighbor graph used ten neighbors; UMAP used
`min_dist=0.5`; Leiden used resolution 0.5 and seed 42. Marker genes were ranked
with Wilcoxon rank-sum tests on the full-gene log-normalized `.raw` matrix and
Benjamini–Hochberg adjustment.

These are descriptive cell-level tests from one donor, not replicated
donor-level differential-expression inference.

## Measured PBMC 3k results

| Measure | Result |
| --- | ---: |
| Cells before / after QC | 2,700 / 2,638 |
| Cells removed | 62 |
| Genes before / after QC | 32,738 / 13,656 |
| Highly variable genes | 2,000 |
| Principal components | 40 |
| Leiden clusters | 6 |

| Cluster | Cells | Fraction |
| --- | ---: | ---: |
| 0 | 1,190 | 45.11% |
| 1 | 341 | 12.93% |
| 2 | 640 | 24.26% |
| 3 | 12 | 0.45% |
| 4 | 419 | 15.88% |
| 5 | 36 | 1.36% |

Clusters 3 and 5 are small and require particular caution.

## CellTypist annotation

The workflow used the fixed `Healthy_COVID19_PBMC.pkl` v1 model and verified
SHA-256
`1eee68d23b4ac0f6999a0a92d6c950d671f31fcdb003ce0760c580b0bd6a3b79`.
The model declares 51 cell types and 3,443 features.

An initial certificate error remains visible in
`logs/annotation_failure.log`. The download path was changed to HTTPS Requests
with its certificate trust bundle; the subsequent checksum and annotation run
succeeded. This retained log is historical, not an unresolved stage failure.

The query contained 2,638 cells and 13,656 genes. Inverse-log totals ranged from
9,999.9989 to 10,000.0015, matching the required normalization. It overlapped
2,818/3,443 model features (81.85%). CellTypist returned 2,638 per-cell
best-match predictions with majority voting disabled and no fallback labels.

The table below is cluster consensus, not annotation accuracy. Mean confidence
is for cells assigned the consensus label, not all cells in the cluster.

| Cluster | Most common prediction | Cells | Fraction | Mean confidence |
| --- | --- | ---: | ---: | ---: |
| 0 | CD4.CM | 549 | 46.13% | 0.154 |
| 1 | B_naive | 270 | 79.18% | 0.559 |
| 2 | CD14_mono | 399 | 62.34% | 0.627 |
| 3 | Platelets | 12 | 100.00% | 1.000 |
| 4 | NK_16hi | 141 | 33.65% | 0.874 |
| 5 | DC2 | 17 | 47.22% | 0.624 |

Weak consensus in clusters 0, 4, and 5 must not be hidden by assigning one
unqualified cluster label.

## Evidence-bounded interpretation

The saved result used the deterministic backend. Metadata records
`validated: true` and `external_api_called: false`. Labels must come from each
cluster's allowlist and every evidence ID must be present in the saved packet.
Cluster 0 remained unresolved; clusters 1–3 received moderate-confidence
CellTypist hypotheses; clusters 4–5 received low-confidence marker-panel
hypotheses with explicit uncertainty.

Passing this evidence contract shows traceability, not biological truth.

## Independent five-class benchmark

### Reference design

The benchmark used official 10x FACS-enriched, Cell Ranger-filtered panels:

| Reference class | Source population | Reported cells | Reported purity | Sampled |
| --- | --- | ---: | ---: | ---: |
| B cell | [CD19+ B cells](https://www.10xgenomics.com/datasets/cd-19-plus-b-cells-1-standard-1-1-0) | 10,085 | 100% | 500 |
| Monocyte | [CD14+ monocytes](https://www.10xgenomics.com/datasets/cd-14-plus-monocytes-1-standard-1-1-0) | 2,612 | 98% | 500 |
| CD4+ T cell | [CD4+ helper T cells](https://www.10xgenomics.com/datasets/cd-4-plus-helper-t-cells-1-standard-1-1-0) | 11,213 | 99% | 500 |
| CD8+ T cell | [CD8+ cytotoxic T cells](https://www.10xgenomics.com/datasets/cd-8-plus-cytotoxic-t-cells-1-standard-1-1-0) | 10,209 | 98% | 500 |
| NK cell | [CD56+ natural killer cells](https://www.10xgenomics.com/datasets/cd-56-plus-natural-killer-cells-1-standard-1-1-0) | 8,385 | 92% | 500 |

All source archives have project-frozen SHA-256 digests that matched the
retrieved bytes. Every cached extracted matrix file was also compared with its
corresponding verified tar member before loading. Sampling was
deterministic with seed 42. The 2,500-cell cohort shared 32,738 genes. Reference
labels were kept in a separate pandas Series and absent from the AnnData object
passed to CellTypist, reducing direct label leakage.

The exact CellTypist-to-broad-class mapping used for scoring is saved with the
run and covered by the artifact manifest. Predictions outside or explicitly
excluded from the mapping were recorded as
`__unmapped__` and counted as errors; metrics were not conditioned on covered
cells.

### Measured metrics

| Metric | Result |
| --- | ---: |
| Cells | 2,500 |
| Mapped cells | 2,414 |
| Unmapped cells | 86 |
| Coverage | 0.9656 |
| Strict accuracy | 0.9328 |
| Strict balanced accuracy | 0.9328 |
| Strict macro-F1 | 0.9478 |
| Strict MCC | 0.9184 |
| CellTypist inference runtime | 1.0578 s |

| Class | Precision | Recall | F1 | Support |
| --- | ---: | ---: | ---: | ---: |
| B cell | 0.9804 | 1.0000 | 0.9901 | 500 |
| Monocyte | 0.9915 | 0.9380 | 0.9640 | 500 |
| CD4+ T cell | 0.8908 | 0.9460 | 0.9176 | 500 |
| CD8+ T cell | 0.9776 | 0.7860 | 0.8714 | 500 |
| NK cell | 0.9980 | 0.9940 | 0.9960 | 500 |

The largest visible error mode was CD8-to-CD4 confusion: 55 of 500 CD8
reference cells mapped to CD4, and 46 CD8 cells were unmapped. Runtime times the
annotator call only; it excludes metric calculation, downloads, cohort
preparation, and PBMC 3k preprocessing.

### Benchmark interpretation boundary

These are real measured results, but the reference is not perfect ground truth.
Panel purity is 92–100%; each biological class is a separate library and may be
a separate donor; class is therefore confounded with library/donor and technical
effects. Cell-level metrics cannot estimate donor-to-donor generalization. The
task covers only five broad blood lineages, not fine subtypes or robust
out-of-reference rejection. The model publication differs from the panel
sources, reducing known direct study overlap but not proving the absence of all
indirect training relationships.

## Tests and artifact integrity

A fresh local test run completed with 58 passed tests in 12.56 seconds. It
emitted two non-fatal warnings: sparse zero-centering in a small synthetic
pipeline test and joblib falling back to logical CPU count.

The post-benchmark manifest records nine inputs and 37 outputs. A fresh
verification rehashed all nine inputs and all 37 outputs and found zero missing
files, size mismatches, checksum mismatches, or path issues. The run's
`execution_history.json` records both configs, their hashes, and the completed
`all` and `benchmark` commands.

## Dashboard status

The Streamlit app reads saved artifacts and includes overview, QC, clusters and
markers, CellTypist, interpretation, benchmark, and provenance pages. It does
not recompute analysis or invent absent metrics. A real local process bound to
`127.0.0.1:8501`; both the Streamlit health endpoint and root HTML responded
successfully before the process was stopped normally. This verifies a local
launch, not deployment hardening or service reliability.

## Reproducibility strengths and gaps

Strengths include fixed input/model artifacts, pinned direct dependencies,
seeded stochastic steps, preserved raw counts, resolved configuration, source
manifests, an external reference-label container, frozen mapping, and output
hashes.

Open gaps remain:

- generated data/results are Git-ignored and need regeneration or archival;
- the exact environment lock is macOS/Python 3.12-specific; no container or
  cross-platform reproduction is recorded;
- Python 3.14 and non-macOS environments are untested;
- no independent final scientific audit has examined the completed benchmark.

## Completion assessment

| Deliverable | Evidence-backed status |
| --- | --- |
| PBMC 3k acquisition and analysis | Complete |
| CellTypist annotation | Complete |
| Deterministic evidence-bounded interpretation | Complete |
| Optional OpenAI implementation | Present; no real API success claimed |
| Streamlit artifact viewer | Implemented; local launch and health check passed |
| Automated tests | 58 passed |
| Five-class FACS benchmark | Complete, with stated scope limits |
| Final independent scientific audit | Not performed |
| Production readiness | Not demonstrated |

See `SCIENTIFIC_LIMITATIONS.md` before using these results in a scientific or
career claim.
