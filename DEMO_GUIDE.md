# CellPilot AI Demo Guide

This guide presents only results verified on disk. It supports a short technical
interview walkthrough and local reproduction.

## Pre-demo checks

```bash
.venv/bin/python --version
.venv/bin/python -m pip check
.venv/bin/cellpilot verify --run-dir results/pbmc3k/verified_run
```

The verified snapshot used Python 3.12.5 and returned `valid: true`, 9 inputs
and 37 outputs checked, and zero issues. Do not use Python 3.14 unless the full
stack has been separately validated there.

## Reproduce the workflow

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install '.[ai,test]'
.venv/bin/python -m pip check
.venv/bin/cellpilot fetch --config configs/pbmc3k.yaml
.venv/bin/cellpilot all --config configs/pbmc3k.yaml
```

For the exact verified macOS/Python 3.12 package set, install
`requirements-lock-py312-macos.txt` followed by
`.venv/bin/python -m pip install --no-deps .` instead of resolving the extras
again.

The full command creates a timestamped run and updates
`results/LATEST_RUN.txt`. To add the real five-class benchmark to a completed
run:

```bash
.venv/bin/cellpilot benchmark \
  --config configs/benchmark.yaml \
  --run-dir results/pbmc3k/verified_run
```

## Launch the dashboard

```bash
CELLPILOT_RUN_DIR=results/pbmc3k/verified_run \
.venv/bin/streamlit run src/cellpilot/dashboard.py
```

Open the local URL Streamlit prints. If it is not reachable, preserve the
terminal output and do not claim a successful live-app demo.

## Seven-minute walkthrough

### 1. Provenance

Open **Provenance & limitations**.

Say:

> The workflow starts from the official 10x PBMC 3k filtered UMI matrix. The
> input archive, the frozen CellTypist model, and five benchmark archives are
> recorded as inputs, while a manifest hashes every generated artifact.

Show the 9-input/37-output verification and the execution history for the
analysis and benchmark configurations.

### 2. QC and normalization

Open **QC**.

Say:

> The analysis began with 2,700 cells and 32,738 genes. Prespecified PBMC 3k QC
> retained 2,638 cells and 13,656 genes. Raw counts are preserved separately;
> normalization targets 10,000 counts per cell and applies log1p once.

Do not call the QC cutoffs universal. Mention that the matrix was already Cell
Ranger filtered and that the workflow does not model doublets or ambient RNA.

### 3. Clustering and markers

Open **Clusters & markers**.

Explain that 2,000 HVGs, 40 PCs, ten neighbors, and Leiden resolution 0.5
produced six clusters with seed 42. The real marker table uses Wilcoxon rank-sum
on full-gene log-normalized data with BH correction.

Say clearly that one donor provides no donor-level biological replication and
that UMAP geometry is not proof of identity.

### 4. CellTypist predictions

Open **CellTypist annotation**.

Say:

> The primary output is per-cell best match with majority voting disabled. The
> fixed model checksum matched, 2,818 of 3,443 model features overlapped, and the
> workflow produced 2,638 predictions without fallback labels.

Use the consensus table to show uncertainty: cluster 1 is 79.18% `B_naive`,
cluster 2 is 62.34% `CD14_mono`, but clusters 0, 4, and 5 have only 46.13%,
33.65%, and 47.22% consensus. These are within-cluster summaries, not accuracy.

### 5. Evidence-grounded interpretation

Open **Evidence-grounded interpretation**.

Say:

> The saved run uses the deterministic interpreter. Each proposed label must be
> in a cluster allowlist, cited evidence IDs must exist, and ambiguous clusters
> can remain unresolved. No external API was called.

Show that cluster 0 is unresolved and that low-confidence marker hypotheses
carry explicit uncertainty. Local schema validation is not biological
validation.

### 6. Independent benchmark

Open **Benchmark**.

Say:

> The benchmark uses 500 cells from each of five official 10x FACS-enriched
> panels. Reference labels are absent from the AnnData passed to CellTypist, the
> exact broad-class mapping used for scoring is saved and hash-covered, and all
> unmapped predictions count as errors.

Present the measured results:

- 2,500 cells and 96.56% mapping coverage;
- 93.28% strict accuracy and balanced accuracy;
- 94.78% strict macro-F1;
- 0.9184 MCC;
- 1.0578 seconds for the CellTypist inference call, excluding metric calculation,
  download, and cohort preparation.

Then show the limitation rather than burying it: each class is a separate
library with reported purity of 92–100%. The metrics do not estimate donor-level
generalization. CD8 was the hardest class (F1 0.8714; recall 0.7860), including
55 CD8-to-CD4 calls and 46 unmapped CD8 cells.

### 7. Close with reproducibility

Distinguish raw data, deterministic computed outputs, model predictions, and
interpretation. The strongest story is not “95% accurate AI”; it is that the
project measures a prespecified task, records every caveat, fails closed on
changed artifacts or malformed inputs, and refuses to turn unmeasured claims
into resume numbers.

## Run tests during a technical demo

```bash
MPLCONFIGDIR=work/cache/matplotlib \
NUMBA_CACHE_DIR=work/cache/numba \
.venv/bin/python -m pytest -q
```

The documented snapshot passed 58 tests in 12.56 seconds with two non-fatal
environment warnings. Report the live result if it differs.

## Optional OpenAI path

The API backend should not be enabled without a valid credential and explicit
intent to make an external request.

For an authorized live test:

1. Set `OPENAI_API_KEY` in the environment; never put it in Git.
2. Copy the configuration and set `interpretation.backend: openai`.
3. Select an explicit supported model.
4. Run `cellpilot interpret` against a completed run.
5. Confirm actual response metadata and local evidence validation.

A mocked integration test validates the software contract only. It is not
evidence that an external API call or scientific interpretation succeeded.

## Static fallback

If Streamlit cannot launch, show the saved QC, UMAP, marker, CellTypist, benchmark
confusion-matrix, interpretation, and manifest artifacts directly. State that
this is a static artifact review, not a successful live application demo.

## Claims to avoid

- Do not say the benchmark provides perfect ground truth.
- Do not generalize its five broad classes to fine subtypes.
- Do not claim donor-level validation or clinical readiness.
- Do not say an LLM generated the verified interpretation.
- Do not say the final independent scientific audit is complete.
