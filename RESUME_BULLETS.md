# CellPilot AI Resume Package

Use only claims supported by the current repository. Keep the benchmark scope
and limitations when quoting its numbers.

## Evidence-safe project title

**CellPilot AI — Reproducible scRNA-seq analysis and annotation workflow**

## Resume bullets supported now

- Built a reproducible Python/Scanpy workflow that checksum-verified and
  processed the public 10x PBMC 3k matrix from 2,700 cells through QC,
  normalization, PCA, UMAP, six-cluster Leiden analysis, and marker ranking,
  retaining 2,638 cells and 13,656 genes.

- Integrated a frozen, SHA-256-verified CellTypist PBMC model to generate 2,638
  per-cell best-match annotations with 81.85% model-feature overlap, strict
  normalization checks, majority voting disabled, and no fabricated fallback
  labels.

- Benchmarked broad-lineage annotation on 2,500 cells sampled from five public
  10x FACS-enriched panels, measuring 0.9328 strict accuracy, 0.9478 macro-F1,
  0.9328 balanced accuracy, 0.9184 MCC, and 0.9656 coverage while counting
  unmapped predictions as errors.

- Designed an evidence-bounded interpretation layer that restricts labels to
  cluster-specific allowlists, validates cited evidence identifiers, and
  preserves unresolved or low-confidence conclusions; the verified run used a
  deterministic backend with no external API call.

- Added reproducibility controls including pinned direct dependencies, seeded
  stochastic analysis, separated raw/normalized matrices, source and model
  checksums, resolved configuration, and a manifest that rehashed nine inputs
  and 37 generated outputs with zero issues.

- Implemented a read-only Streamlit interface for real QC, UMAP, markers,
  annotation, interpretation, benchmark, and provenance artifacts, with 58
  passing local tests covering scientific invariants, failure behavior,
  manifest integrity, and benchmark label isolation.

## Required benchmark qualifier

Use this sentence in a portfolio README, LinkedIn post, or interview when space
allows:

> The benchmark is a five-class descriptive evaluation on separate 10x
> FACS-enriched libraries with reported purity of 92–100%; it does not establish
> donor-level generalization, fine-subtype accuracy, or clinical validity.

The benchmark CellTypist inference call took 1.0578 seconds, but that excludes
metric calculation, cohort preparation, data download, and PBMC 3k processing.
Do not use it as an end-to-end speed claim.

## Short LinkedIn description

CellPilot AI is an alpha computational-biology project that turns public
single-cell count matrices into traceable Scanpy analysis, CellTypist annotation,
evidence-linked interpretation, and a measured broad-lineage benchmark. The
verified PBMC 3k run processes 2,700 cells, preserves raw and normalized
representations separately, records model and data checksums, and reports 0.9328
strict accuracy on a 2,500-cell five-class FACS-enriched benchmark while
explicitly documenting purity and separate-library confounding.

## Interview summary

> I built CellPilot AI to demonstrate that an AI-assisted biology workflow can
> be reproducible and honest about uncertainty. The deterministic pipeline
> starts from a checksum-verified public count matrix, preserves raw counts,
> applies normalization once, and hashes its outputs. Annotation uses a fixed
> CellTypist model with input and feature-overlap checks. The interpretation
> layer cannot cite evidence or propose labels outside the saved packet. For the
> benchmark, I kept reference labels outside the model input, froze the class
> mapping, and counted unmapped outputs as errors. I also explain why the
> measured 93.28% strict accuracy does not imply donor-level or clinical
> validity.

## Strong technical discussion points

- Why full-gene log-normalized input is kept separate from scaled PCA input.
- Why cluster consensus cannot substitute for independent reference labels.
- Why direct label leakage is reduced by keeping references outside AnnData.
- Why strict metrics include unmapped predictions instead of conditioning on
  coverage.
- Why CD8 performance (F1 0.8714, recall 0.7860) is more informative than only
  quoting overall accuracy.
- Why FACS purity and one-library-per-class confounding constrain the benchmark.
- Why schema/evidence validation prevents some hallucinations but cannot prove
  biological truth.

## Claims not supported

Do not write or say:

- “95% accurate cell-type annotation” without naming the five-class FACS task
  and strict metric; macro-F1 is not accuracy.
- “Validated across donors,” “independent biological replication,” or
  “generalizes to disease cohorts.”
- “Fine-grained cell-subtype accuracy” or “out-of-reference detection.”
- “Reduced analysis time by X%” or “saved X hours.”
- “Production-ready,” “clinically validated,” or “deployment-ready.”
- “Used GPT/OpenAI for the verified interpretations.” The saved result records
  `external_api_called: false`.
- “Passed a final independent scientific audit.” That audit has not occurred.

## Suggested skills line

Python 3.12, Scanpy, AnnData, CellTypist, scikit-learn, Streamlit, Pydantic,
pytest, scRNA-seq QC, dimensionality reduction, clustering, marker analysis,
annotation benchmarking, model provenance, scientific reproducibility

List OpenAI API execution only after a real authorized request succeeds and its
validated metadata is saved; mocked integration is not equivalent.
