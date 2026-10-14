# CellPilot AI — Detailed Usage Guide

> This is the original technical README from verified commit `2c7ba636`, kept
> verbatim below this note. Paths are relative to the repository root. The
> portfolio overview is in the top-level [README](../README.md).

CellPilot AI is a reproducible single-cell RNA-seq project built on real public
data. It acquires the 10x Genomics PBMC 3k count matrix, verifies inputs,
performs deterministic Scanpy analysis, annotates cells with a frozen CellTypist
model, creates evidence-bounded biological hypotheses, benchmarks broad-lineage
annotation, and serves saved artifacts in a read-only Streamlit dashboard.

This is an alpha research project, not a clinical tool or a production-ready
analysis platform. Labels and interpretations remain hypotheses requiring
biological review.

## Verified results

The following are backed by files in `results/pbmc3k/verified_run` as of
2026-10-14:

| Item | Measured result |
| --- | --- |
| PBMC 3k input | 2,700 cells × 32,738 genes; configured archive SHA-256 matched |
| QC output | 2,638 cells and 13,656 genes retained; 62 cells removed |
| Analysis | 2,000 HVGs, 40 PCs, six Leiden clusters, seed 42 |
| Normalization | `normalize_total(target_sum=10000)` then `log1p` exactly once |
| CellTypist | 2,638 best-match predictions; 2,818/3,443 model features overlapped (81.85%) |
| Interpretation | Deterministic backend validated against saved evidence; no external API call |
| FACS benchmark | 2,500 cells, five broad classes; 0.9328 strict accuracy; 0.9478 macro-F1; 0.9184 MCC |
| Artifact integrity | Manifest verified 9 inputs and 37 outputs with zero issues |
| Local tests | 58 passed in 12.56 seconds; two non-fatal environment warnings |
| Final scientific audit | **Not performed:** the existing audit predates this implementation |

Benchmark results are descriptive cell-level measurements on five separate 10x
FACS-enriched libraries with reported purity of 92–100%. They do not demonstrate
donor-level generalization or fine-subtype accuracy. Read
`SCIENTIFIC_LIMITATIONS.md` before interpreting them.

## Project layout

```text
configs/                         Analysis and benchmark settings
src/cellpilot/                   Pipeline, annotation, benchmark, interpretation, dashboard
src/cellpilot/resources/         Marker evidence and frozen benchmark label mapping
tests/                           Local unit and integration-style tests
data/                            Downloaded data/models; ignored by Git
results/                         Generated scientific artifacts; ignored by Git
PROJECT_REPORT.md                Verified methods and measured results
DEMO_GUIDE.md                    Reproduction and demo walkthrough
SCIENTIFIC_LIMITATIONS.md        Scope, failure modes, and untested claims
RESUME_BULLETS.md                Evidence-safe project wording
requirements-lock-py312-macos.txt Exact verified macOS/Python 3.12 environment
```

## Installation

Python 3.12 is required. The verified run used CPython 3.12.5; the host's Python
3.14 interpreter was intentionally not used because this scientific stack was
not validated there.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install '.[ai,test]'
.venv/bin/python -m pip check
```

Use a normal local install, as above. An editable install was unreliable on this
host because its Python site machinery skipped the hidden editable-install
`.pth` file. Direct dependencies are pinned in `pyproject.toml`. The verified
macOS/Python 3.12 environment is captured in
`requirements-lock-py312-macos.txt`; it is not presented as a cross-platform
lock.

To recreate that exact platform-specific package set instead of resolving
transitive dependencies again:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock-py312-macos.txt
.venv/bin/python -m pip install --no-deps .
```

## Run the PBMC 3k workflow

The deterministic workflow requires no API key.

```bash
.venv/bin/cellpilot fetch --config configs/pbmc3k.yaml
.venv/bin/cellpilot all --config configs/pbmc3k.yaml
```

The full command creates a timestamped run and updates
`results/LATEST_RUN.txt`. Each completed run contains a processed AnnData object,
real plots and tables, CellTypist probabilities, interpretation evidence, stage
status, and an integrity manifest.

Verify the saved run without rerunning analysis:

```bash
.venv/bin/cellpilot verify --run-dir results/pbmc3k/verified_run
```

The current snapshot returned `valid: true`, rehashed 9 recorded inputs and 37
outputs, and reported no integrity issues. `execution_history.json` separately
records the completed `all` and `benchmark` invocations and their config hashes.

## Run the benchmark

```bash
.venv/bin/cellpilot benchmark \
  --config configs/benchmark.yaml \
  --run-dir results/pbmc3k/verified_run
```

The completed run deterministically sampled 500 cells from each of five official
10x FACS-enriched PBMC panels. Reference labels were stored outside the AnnData
object given to CellTypist. The exact broad-class mapping used for scoring was
saved with the run and hash-covered; unmapped predictions counted as errors in
every strict metric.

| Metric | Value |
| --- | ---: |
| Cells | 2,500 |
| Mapped / coverage | 2,414 / 0.9656 |
| Strict accuracy | 0.9328 |
| Strict balanced accuracy | 0.9328 |
| Strict macro-F1 | 0.9478 |
| Strict MCC | 0.9184 |
| CellTypist inference runtime | 1.0578 s |

Runtime times the CellTypist annotator call only, not metric calculation, dataset
download, cohort preparation, or the PBMC 3k analysis. See `PROJECT_REPORT.md`
for class-level results and limitations.

## Run tests

```bash
MPLCONFIGDIR=work/cache/matplotlib \
NUMBA_CACHE_DIR=work/cache/numba \
.venv/bin/python -m pytest -q
```

The documented execution passed 58 tests. They cover count validation,
normalization and count-layer preservation, CellTypist input/model checks,
evidence validation, manifest tamper detection, and benchmark label isolation.
They do not replace independent biological review.

## Launch the dashboard

```bash
CELLPILOT_RUN_DIR=results/pbmc3k/verified_run \
.venv/bin/streamlit run src/cellpilot/dashboard.py
```

The dashboard displays overview, QC, clusters and markers, annotation,
interpretation, benchmark, and provenance pages. It reads saved artifacts and
does not rerun the scientific pipeline. A real local launch on
`127.0.0.1:8501` returned `ok` from Streamlit's health endpoint; this is a local
smoke test, not a public deployment claim.

## Optional OpenAI interpretation

The verified interpretation used the deterministic backend and records
`external_api_called: false`. The optional backend uses Responses API structured
output and then rejects unknown evidence IDs, unsupported labels, duplicate
clusters, and missing clusters locally.

To use it, set `interpretation.backend: openai` in a configuration file and
provide `OPENAI_API_KEY` through the environment. Never commit a key. A mocked
integration test is not evidence that a real API request succeeded.

## Primary sources

- [10x Genomics PBMC 3k dataset](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-0-0)
- [Scanpy documentation](https://scanpy.readthedocs.io/)
- [CellTypist documentation](https://celltypist.readthedocs.io/)
- [CellTypist model source publication](https://doi.org/10.1038/s41591-021-01329-2)
- [OpenAI Structured Outputs](https://developers.openai.com/api/docs/guides/structured-outputs?api-mode=responses)

Never treat generated prose, clusters, marker rankings, model predictions, or
benchmark agreement as a substitute for experimental validation.
