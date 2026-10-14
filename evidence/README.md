# Evidence bundle

Small, reviewable artifacts copied from the verified run
`results/pbmc3k/verified_run` at commit `2c7ba636`. The full run directory
(processed AnnData, per-cell probability tables, raw downloads, and the
CellTypist model) is intentionally **not** committed; regenerate it with the
commands in the [README](../README.md#reproduce-the-results).

| Folder | Contents |
| --- | --- |
| `analysis/` | QC/clustering summary, cluster sizes, CellTypist cluster consensus and model metadata, evidence-bounded interpretation (JSON + Markdown) |
| `benchmark/` | Strict metrics, per-class metrics, raw and row-normalized confusion matrices, frozen 51-label → 5-class mapping, FACS data manifest with SHA-256 checksums |
| `provenance/` | Run integrity manifest (9 inputs / 37 outputs), execution history, PBMC 3k source manifest, stage status |

Figures from the same run are in [`docs/figures/`](../docs/figures/).

## Redaction note

Two files recorded absolute paths from the original development machine. For
publication, the local project-root prefix was replaced with the literal token
`<PROJECT_ROOT>`; no other bytes were changed:

- `provenance/manifest.json`
- `analysis/celltypist_metadata.json`

Because of this redaction, these two files will not match the hashes recorded
for them inside the original run. SHA-256 values of the **unredacted** originals
(from `cellpilot-ai-evidence-2c7ba63.zip`, SHA-256
`65a138457f293f298823eb55c6ec7971b712587118ff95d9c72d1a4df1633f72`) are listed
below so the provenance chain stays auditable. All other files here are
byte-identical to the evidence bundle.

```text
8fa27332066c208b8d411625c208e6eab5b4e97e7b5d9bc9259164517f357f07  analysis/celltypist_metadata.json
d36f58bb7e78655bb3ac1f064e237a0be337a872b040981be0f6ec7f4773c73d  provenance/manifest.json
```
