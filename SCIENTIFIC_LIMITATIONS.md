# Scientific Limitations and Claim Boundaries

This document identifies what CellPilot AI has and has not demonstrated as of
the 2026-10-14 verified run.

## Intended use

CellPilot AI is an educational and portfolio-scale exploratory research
workflow. It is not intended for diagnosis, treatment decisions, clinical
reporting, regulated use, or unsupervised publication-quality inference.

## Input data

- PBMC 3k represents one healthy donor. It cannot establish reproducibility
  across donors, demographics, diseases, sites, batches, or protocols.
- The source is a Cell Ranger-filtered gene-barcode matrix, not raw droplets or
  FASTQ. Alignment, barcode calling, and upstream filtering are outside scope.
- A filtered matrix cannot support complete empty-droplet or ambient-RNA
  analysis. The workflow does not explicitly call doublets or model ambient RNA.
- PBMC findings should not be generalized to other tissues or species.

## Quality control

- The 200–2,500 detected-gene cutoff, less than 5% mitochondrial threshold, and
  three-cell gene filter were prespecified for this demonstration. They are not
  universal best practices.
- No threshold-sensitivity analysis was performed. The six-cluster result may
  change under other QC choices.
- Total-count cutoffs were not applied, and no dedicated doublet, stress, or
  library-complexity assessment was performed.
- Removing only 62 cells is not itself evidence that QC was biologically
  optimal.

## Normalization and embeddings

- Library-size normalization to 10,000 followed by `log1p` is an exploratory
  choice, not a universally correct UMI model.
- Alternatives such as Pearson residuals or model-based normalization were not
  compared.
- UMAP is a neighbor-graph visualization. Distances, gaps, and cluster shapes
  must not be interpreted as quantitative biological evidence.
- A fixed seed improves repeatability but does not establish robustness. No
  multi-seed or parameter-stability study is reported.

## Clustering and marker statistics

- Leiden resolution 0.5 yielding six clusters is one parameterization, not a
  known true number of populations.
- Clusters 3 and 5 contain only 12 and 36 cells. Their markers and labels are
  especially sensitive to noise and sampling.
- Marker tests use cells from one donor as observational units. They do not
  provide donor-level replication and can create pseudoreplication if described
  as population-level differential expression.
- BH adjustment addresses multiplicity under its assumptions; it cannot repair
  absent biological replication or circular selection from clustering the same
  data.

## CellTypist annotation

- CellTypist predictions are model outputs, not reference labels.
- An immune-cell model is domain-relevant but not guaranteed valid for PBMC 3k.
- The query overlapped 2,818 of 3,443 model features (81.85%); missing features
  may affect predictions.
- Majority voting was disabled. Cluster consensus is descriptive only.
- Consensus was weak in clusters 0 (46.13%), 4 (33.65%), and 5 (47.22%). High
  confidence among cells assigned one label does not make it representative of
  the whole cluster.
- No fallback labels are created. This avoids fabrication but deliberately
  exposes failure or ambiguity.
- CellTypist confidence calibration has not been evaluated.
- Training-data overlap and indirect study leakage have not been exhaustively
  excluded.

## Interpretation

- The verified interpretation is deterministic and rule-based. No external LLM
  was called.
- Restricting labels and evidence IDs prevents some hallucinations but cannot
  make incomplete evidence complete or guarantee a correct label.
- Marker panels are simplified broad-lineage heuristics and cannot represent
  every subtype, activation state, or artifact.
- The interpreter does not independently retrieve and verify a primary
  literature citation for each cluster conclusion.
- `validated: true` means a local schema/evidence contract passed, not that the
  biology was validated.
- The optional OpenAI path uses structured output plus local checks. Mocked API
  tests validate integration only, not a real request or scientific quality.

## Five-class benchmark

The benchmark produced real descriptive metrics on 2,500 cells: 0.9328 strict
accuracy and balanced accuracy, 0.9478 macro-F1, 0.9184 MCC, and 0.9656 mapping
coverage. These values have important boundaries:

- Reference labels are FACS-enriched population names, not perfect per-cell
  ground truth. Source pages report 92–100% purity.
- Each class is a separate 10x library and may be a separate donor. Biological
  class is confounded with library/donor and technical effects.
- Cells within a panel are not independent biological replicates. Cell-level
  metrics do not estimate donor-to-donor generalization.
- The task covers only five broad lineages: B cells, monocytes, CD4 T cells, CD8
  T cells, and NK cells. It does not measure fine-subtype performance.
- A frozen mapping collapses detailed CellTypist predictions. Out-of-scope
  labels count as errors, but broad collapsing can still hide subtype errors.
- Sampling 500 cells per class balances the task but discards data and makes the
  exact result conditional on seed 42.
- Bootstrap confidence intervals were disabled because resampling cells would
  not supply missing biological replication.
- Direct label leakage was reduced by keeping reference labels outside AnnData,
  but this does not rule out all indirect training relationships.
- The model publication differs from these 10x panels, reducing known direct
  study overlap but not proving complete independence.
- CD8 was the weakest class: recall 0.7860 and F1 0.8714, with 55 CD8-to-CD4
  calls and 46 unmapped CD8 cells.
- The 1.0578-second runtime times only the CellTypist annotator call. It excludes
  metric calculation, cohort preparation, downloads, and the primary PBMC
  pipeline and is not an end-to-end performance claim.

## Reproducibility

- Fixed checksums detect changed saved inputs, but remote hosts may disappear.
- A valid manifest verifies bytes, not scientific correctness.
- The exact environment snapshot is locked only for macOS with Python 3.12; no
  cross-platform lock or container is present.
- The verified platform is arm64 macOS. Linux, Windows, x86_64, and Python 3.14
  remain untested.
- Numerical libraries can produce platform-dependent differences despite a
  fixed seed.
- Data, models, and results are Git-ignored and need separate archival or
  regeneration.

## Software and dashboard

- Fifty-eight passing tests do not prove the absence of defects.
- Synthetic and component tests do not replace biological validation.
- The Streamlit app is a local read-only viewer and is not hardened for
  untrusted multi-user deployment.
- A visible page does not prove a pipeline stage succeeded; stage status and the
  manifest must also be checked.
- No service reliability, load, accessibility, browser compatibility, or public
  deployment claim has been established.

## Audit status

`outputs/SCIENTIFIC_AUDIT.md` predates the implementation. It documented that
there was then nothing to audit; it is not a final audit of the current code,
PBMC result, or benchmark.

The requested final independent scientific audit has not been performed. It
should independently examine normalization, marker statistics, label mapping,
data leakage, annotation errors, unsupported interpretation, reproducibility,
and benchmark claims.

## Safe conclusion

Supported:

> CellPilot AI is a functioning alpha workflow that processes one public PBMC
> dataset, records real outputs, applies a fixed annotation model, constrains
> interpretation to saved evidence, and reports a measured five-class
> FACS-enriched benchmark with explicit limitations.

Not supported:

> CellPilot AI is an accurate across-donor, generalizable, production-ready, or
> clinically validated cell-annotation system.
