# CellPilot AI Build, Run, and Test Report

> Historical build report for verified commit `2c7ba636`, reproduced unchanged
> below. Item 7 under *Remaining blockers* (no GitHub remote) has since been
> resolved by publishing this repository; all other statements stand as written.

Report date: 2026-10-14  
Project state: working alpha research workflow  
Git commit: `2c7ba636d50218fdf378b94eb6ac87ce48e6b367` (`main`, clean)  
Evidence run: `results/pbmc3k/verified_run`

## Outcome

CellPilot AI was continued from the existing repository and completed through
Build → Run → Test. The existing PBMC 3k outputs were preserved. A checksum-
verified CellTypist model now annotates the processed cells, the end-to-end run
is complete, the Streamlit application launches against the real artifacts, an
independent five-class annotation benchmark has executed, and the deterministic
and optional AI interpretation paths have been validated within their stated
limits.

This report does **not** claim production readiness, clinical validity, or a
completed independent scientific audit. The final audit was intentionally not
performed in this phase.

## Verified environment

- CPython 3.12.5 in `.venv`; the host Python 3.14.4 was not used.
- Scanpy 1.12.4, AnnData 0.13.4, CellTypist 1.7.1, Streamlit 1.65.0.
- NumPy 2.5.4, pandas 3.0.6, SciPy 1.18.1, scikit-learn 1.7.2.
- Pydantic 2.14.1 and OpenAI Python library 2.54.0.
- `pip check`: `No broken requirements found.`
- Exact platform lock: `requirements-lock-py312-macos.txt`.

## End-to-end PBMC 3k run

The source is the official 10x Genomics PBMC 3k Cell Ranger-filtered count
matrix. The downloaded archive matched the frozen SHA-256
`847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5`.

| Measurement | Actual result |
| --- | ---: |
| Input matrix | 2,700 cells × 32,738 genes |
| Cells after QC | 2,638 |
| Genes after QC | 13,656 |
| Highly variable genes | 2,000 |
| Principal components | 40 |
| Leiden clusters | 6 |
| Random seed | 42 |

Normalization was `normalize_total(target_sum=10000)` followed by `log1p`
exactly once. Filtered integer counts were preserved in `layers["counts"]`, and
marker ranking used Wilcoxon tests on the full-gene log-normalized `.raw`
matrix with Benjamini–Hochberg-adjusted p-values.

## CellTypist resolution and annotation

The initial SSL certificate failure was preserved in the historical run log.
The model downloader was changed to verified HTTPS Requests/certifi handling,
then tested with a new empty download directory. The fresh file matched the
expected 816,569-byte size and SHA-256
`1eee68d23b4ac0f6999a0a92d6c950d671f31fcdb003ce0760c580b0bd6a3b79`.

The final run used `Healthy_COVID19_PBMC.pkl` v1 and produced 2,638 per-cell
best-match predictions. Majority voting was disabled and no fallback labels
were used. The query overlapped 2,818 of the model's 3,443 features (81.85%).
Inverse-log totals were approximately 10,000 for every cell, validating the
normalization expected by the model.

## Independent annotation benchmark

The benchmark used five official 10x FACS-enriched, Cell Ranger-filtered PBMC
panels. It sampled 500 cells per panel with seed 42, for 2,500 cells and 32,738
common genes. Reference labels were kept outside the AnnData object passed to
CellTypist. A frozen mapping accounted for all 51 model labels; excluded or
unmapped predictions counted as errors rather than being dropped.

| Metric | Actual result |
| --- | ---: |
| Cells | 2,500 |
| Mapped / unmapped | 2,414 / 86 |
| Coverage | 0.9656 |
| Strict accuracy | 0.9328 |
| Strict balanced accuracy | 0.9328 |
| Strict macro-F1 | 0.9478 |
| Strict MCC | 0.9184 |
| CellTypist annotator-call runtime | 1.0578 s |

| Reference class | Precision | Recall | F1 | Support |
| --- | ---: | ---: | ---: | ---: |
| B cell | 0.9804 | 1.0000 | 0.9901 | 500 |
| Monocyte | 0.9915 | 0.9380 | 0.9640 | 500 |
| CD4+ T cell | 0.8908 | 0.9460 | 0.9176 | 500 |
| CD8+ T cell | 0.9776 | 0.7860 | 0.8714 | 500 |
| NK cell | 0.9980 | 0.9940 | 0.9960 | 500 |

The largest visible error mode was CD8-to-CD4 confusion: 55 of 500 CD8
reference cells were scored as CD4, and 46 were unmapped. The measured runtime
covers only the CellTypist annotator call, not downloads, data preparation,
metrics, or the main PBMC analysis.

These FACS labels are enriched population labels with reported panel purities
of 92–100%, not perfect cell-level ground truth. Each class is a separate
library and may be a separate donor, so these metrics do not measure donor-level
generalization.

## Interpretation validation

The saved biological interpretation used the deterministic backend and records
`validated: true` and `external_api_called: false`. Six cluster records were
generated with actual cell counts, evidence IDs, constrained labels, confidence
levels, and explicit uncertainty. Cluster 0 remained unresolved rather than
receiving an unsupported label.

The optional OpenAI path uses Responses API structured output and Pydantic
validation, followed by local checks for allowed labels, evidence IDs,
duplicates, and missing clusters. Mocked tests exercised successful parsing and
rejection behavior. No real API key was supplied and no live API success is
claimed.

## Test and application verification

- Final test rerun: **58 passed, 2 warnings in 14.58 seconds**.
- The warnings were non-fatal: sparse zero-centering in a small synthetic test
  and joblib falling back to the logical CPU count.
- Dashboard tests exercised all seven artifact pages.
- A real Streamlit process started on `127.0.0.1:8501`.
- `/_stcore/health` returned `ok`, and `/` returned the Streamlit HTML shell.
- The server then stopped normally.
- Final artifact verification rehashed **9 inputs and 37 outputs**, returning
  `valid: true` with an empty issue list.

## Successfully executed commands

The persisted execution history records:

```text
.venv/bin/cellpilot all --config configs/pbmc3k.yaml --run-dir results/pbmc3k/verified_run
.venv/bin/cellpilot benchmark --config configs/benchmark.yaml --run-dir results/pbmc3k/verified_run
.venv/bin/cellpilot interpret --config configs/pbmc3k.yaml --run-dir results/pbmc3k/verified_run
```

Final verification also executed:

```text
.venv/bin/cellpilot verify --run-dir results/pbmc3k/verified_run
.venv/bin/python -m pytest -q
.venv/bin/python -m pip check
```

The source was installed as a normal local wheel rather than an editable
installation. The built wheel used for the verified installation had SHA-256
`df262852b5b716fee0b41ba2021e240e8285fe677699ca6dc90a15be502de1b0`.

## Repository structure

```text
configs/                         Frozen pipeline and benchmark settings
src/cellpilot/                   Analysis, annotation, benchmark, AI, dashboard
src/cellpilot/resources/         Marker evidence and frozen label mapping
tests/                           Unit and integration-style tests
data/                            Downloaded inputs/models; Git-ignored
results/                         Full generated run; Git-ignored
outputs/                         User-facing report, evidence, plots, source ZIP
README.md                        Install and usage guide
PROJECT_REPORT.md                Methods and measured scientific results
DEMO_GUIDE.md                    Reproduction and interview demo guide
SCIENTIFIC_LIMITATIONS.md        Boundaries and untested claims
RESUME_BULLETS.md                Evidence-safe resume wording
requirements-lock-py312-macos.txt Exact verified platform lock
```

## Packaged deliverables

- `cellpilot-ai-source-2c7ba63.zip`: Git-tracked source at the verified commit.
  SHA-256:
  `2757317e2395aac7142e22d3b99e471152d5f8210f409684fac65b2ae3dec20b`.
- `PROJECT_REPORT.md`: scientific methods, results, and evidence boundaries.
- `DEMO_GUIDE.md`: exact reproduction and demonstration workflow.
- `SCIENTIFIC_LIMITATIONS.md`: known scientific and engineering limitations.
- `RESUME_BULLETS.md`: measured, defensible project statements.
- `demo/`: actual QC, HVG, PCA, UMAP, annotation, and marker plots.
- `benchmark/`: actual metrics, confusion matrices, mapping, and data manifest.
- `analysis/`: PBMC summary, annotation metadata, and interpretation artifacts.
- `provenance/`: run manifest, execution history, source manifest, stage state.

The complete processed AnnData and larger probability/cell tables remain under
`results/pbmc3k/verified_run`; they are intentionally excluded from Git and the
small review bundle.

## Remaining blockers and limitations

1. A final skeptical independent scientific audit has not yet been performed.
2. The optional OpenAI backend has mocked integration coverage but no live API
   execution; no external API result should be claimed.
3. The benchmark covers only five broad PBMC lineages and has library/donor
   confounding, imperfect FACS purity, and no biological replicates per class.
4. The PBMC demonstration uses one donor and filtered matrices; it does not
   model ambient RNA, doublets, batch effects, or condition effects.
5. The exact environment lock is specific to Python 3.12 on macOS. Python 3.14,
   other operating systems, containers, and clean-machine reproduction remain
   untested.
6. The dashboard was validated locally but was not publicly deployed or load-
   tested.
7. No GitHub remote is configured, so the deliverable is a local Git commit and
   source ZIP rather than a pushed GitHub repository.
8. `outputs/SCIENTIFIC_AUDIT.md` is the historical pre-build audit. It correctly
   reported that no implementation existed at that time and must not be treated
   as an audit of this completed build.

## Primary references

- [10x Genomics PBMC 3k](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-0-0)
- [CellTypist documentation](https://celltypist.readthedocs.io/)
- [CellTypist model publication](https://doi.org/10.1038/s41591-021-01329-2)
- [OpenAI Structured Outputs](https://developers.openai.com/api/docs/guides/structured-outputs?api-mode=responses)
- [OpenAI Python API reference](https://developers.openai.com/api/reference/python)

