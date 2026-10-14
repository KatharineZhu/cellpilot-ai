"""Command-line interface for reproducible CellPilot runs."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import pandas as pd
import yaml


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _default_run_id() -> str:
    return datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _read_json(path: Path) -> Any | None:
    """Read a JSON artifact when present, otherwise return ``None``."""

    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _load_status(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "stage_status.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {"run_id": run_dir.name, "started_at_utc": _utc_now(), "stages": {}}


def _set_stage(
    run_dir: Path,
    stage: str,
    status: str,
    *,
    detail: str | None = None,
) -> None:
    payload = _load_status(run_dir)
    record = payload["stages"].setdefault(stage, {})
    record["status"] = status
    record["updated_at_utc"] = _utc_now()
    if detail:
        record["detail"] = detail
    else:
        record.pop("detail", None)
    _write_json(run_dir / "stage_status.json", payload)


def _run_stage(run_dir: Path, stage: str, action: Callable[[], Any]) -> Any:
    _set_stage(run_dir, stage, "running")
    try:
        result = action()
    except Exception as exc:
        _set_stage(run_dir, stage, "failed", detail=f"{type(exc).__name__}: {exc}")
        log_dir = run_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        (log_dir / f"{stage}_failure.log").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
        raise
    _set_stage(run_dir, stage, "completed")
    return result


def _resolve_project_root(value: str) -> Path:
    root = Path(value).expanduser().resolve()
    if not (root / "pyproject.toml").is_file():
        raise FileNotFoundError(f"Not a CellPilot project root: {root}")
    return root


def _resolve_config(root: Path, value: str) -> Path:
    path = Path(value)
    path = path if path.is_absolute() else root / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def _resolve_run_dir(root: Path, configured_output: Path, value: str | None) -> Path:
    if value:
        path = Path(value)
        resolved = (path if path.is_absolute() else root / path).resolve()
    else:
        output = configured_output
        output = output if output.is_absolute() else root / output
        resolved = (output / _default_run_id()).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError(f"Run directory must stay inside project root: {resolved}")
    return resolved


def _write_resolved_config(run_dir: Path, config: Any) -> None:
    payload = config.model_dump(mode="json")
    path = run_dir / "config.resolved.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def _record_execution(
    run_dir: Path,
    root: Path,
    *,
    operation: str,
    config_path: Path,
) -> None:
    """Append one completed CLI operation before the run manifest is rebuilt."""

    from .provenance import sha256_file

    history_path = run_dir / "execution_history.json"
    payload = _read_json(history_path) or {"schema_version": "1.0", "events": []}
    config_resolved = config_path.resolve()
    try:
        config_record = str(config_resolved.relative_to(root.resolve()))
    except ValueError:
        config_record = str(config_resolved)
    payload["events"].append(
        {
            "recorded_at_utc": _utc_now(),
            "operation": operation,
            "status": "completed",
            "command": list(sys.argv),
            "config_path": config_record,
            "config_sha256": sha256_file(config_resolved),
        }
    )
    _write_json(history_path, payload)


def _collect_run_inputs(root: Path, run_dir: Path) -> list[Path]:
    """Collect every config and immutable external artifact used by this run."""

    from .provenance import sha256_file

    candidates: list[Path] = []
    history = _read_json(run_dir / "execution_history.json") or {}
    for event in history.get("events", []):
        value = event.get("config_path") if isinstance(event, dict) else None
        if not value:
            continue
        path = Path(str(value))
        resolved_config = (path if path.is_absolute() else root / path).resolve()
        recorded_sha = event.get("config_sha256")
        if not resolved_config.is_file():
            raise FileNotFoundError(f"recorded run config is missing: {resolved_config}")
        if recorded_sha and sha256_file(resolved_config) != recorded_sha:
            raise ValueError(
                "a configuration used by an earlier run operation has changed: "
                f"{resolved_config}"
            )
        candidates.append(resolved_config)

    if (run_dir / "analysis_summary.json").is_file():
        candidates.append(
            root / "data/raw/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
        )
    if (run_dir / "annotation/celltypist_metadata.json").is_file():
        candidates.append(root / "data/celltypist/Healthy_COVID19_PBMC.pkl")
    if (run_dir / "benchmark/metrics.json").is_file():
        candidates.extend(
            sorted((root / "data/benchmark/facs_5class/archives").glob("*.tar.gz"))
        )

    unique: dict[Path, None] = {}
    for candidate in candidates:
        resolved = candidate.resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"recorded run input is missing: {resolved}")
        unique[resolved] = None
    return sorted(unique, key=lambda path: str(path))


def _refresh_manifest(
    root: Path,
    run_dir: Path,
    *,
    config_path: Path,
    random_seed: int,
) -> None:
    """Create one unified manifest after the latest run operation."""

    from .provenance import write_manifest

    summary = _read_json(run_dir / "analysis_summary.json") or {}
    warnings = list(summary.get("limitations", []))
    if (run_dir / "benchmark/metrics.json").is_file():
        from .benchmark import BENCHMARK_LIMITATIONS

        warnings.extend(BENCHMARK_LIMITATIONS)
    _set_stage(run_dir, "manifest", "running")
    _set_stage(run_dir, "manifest", "completed")
    try:
        write_manifest(
            run_dir,
            project_root=root,
            inputs=_collect_run_inputs(root, run_dir),
            config_path=config_path,
            command=sys.argv,
            random_seed=random_seed,
            warnings=warnings,
        )
    except Exception as exc:
        _set_stage(
            run_dir,
            "manifest",
            "failed",
            detail=f"{type(exc).__name__}: {exc}",
        )
        raise


def _save_annotation_artifacts(
    adata: Any,
    result: Any,
    run_dir: Path,
) -> pd.DataFrame:
    import matplotlib.pyplot as plt
    import scanpy as sc

    predictions = result.predictions.copy()
    predictions.to_csv(run_dir / "tables" / "celltypist_cells.csv")
    result.probabilities.to_csv(
        run_dir / "tables" / "celltypist_probabilities.csv.gz", compression="gzip"
    )
    adata.obs["celltypist_label"] = predictions.loc[
        adata.obs_names.astype(str), "predicted_label"
    ].to_numpy()
    adata.obs["celltypist_confidence"] = predictions.loc[
        adata.obs_names.astype(str), "confidence"
    ].to_numpy(dtype=float)

    grouped = (
        adata.obs.reset_index(names="cell_id")
        .groupby(["leiden", "celltypist_label"], observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            mean_confidence=("celltypist_confidence", "mean"),
        )
        .reset_index()
    )
    totals = grouped.groupby("leiden", observed=True)["n_cells"].transform("sum")
    grouped["fraction"] = grouped["n_cells"] / totals
    consensus = (
        grouped.sort_values(
            ["leiden", "fraction", "mean_confidence"],
            ascending=[True, False, False],
        )
        .groupby("leiden", observed=True)
        .head(1)
        .rename(columns={"leiden": "cluster", "celltypist_label": "predicted_label"})
        .reset_index(drop=True)
    )
    consensus.to_csv(
        run_dir / "tables" / "celltypist_cluster_consensus.csv", index=False
    )

    annotation_dir = run_dir / "annotation"
    annotation_dir.mkdir(parents=True, exist_ok=True)
    _write_json(
        annotation_dir / "celltypist_metadata.json",
        {
            "status": "completed",
            "model": result.model_metadata,
            "input_validation": result.input_validation,
            "n_predictions": int(len(predictions)),
            "primary_prediction": "best match per cell; majority voting disabled",
        },
    )
    sc.pl.umap(
        adata,
        color="celltypist_label",
        frameon=False,
        title="PBMC 3k — CellTypist predictions",
        show=False,
    )
    plt.gcf().savefig(
        run_dir / "plots" / "umap_celltypist.png", dpi=160, bbox_inches="tight"
    )
    plt.close("all")
    adata.write_h5ad(run_dir / "data" / "pbmc3k_processed.h5ad", compression="gzip")
    return consensus


def _run_interpretation(root: Path, config: Any, run_dir: Path) -> None:
    from .interpretation import (
        build_evidence_packet,
        deterministic_interpret,
        load_marker_registry,
        openai_interpret,
        write_interpretation_artifacts,
    )

    markers = pd.read_csv(run_dir / "tables" / "markers.csv")
    counts = pd.read_csv(run_dir / "tables" / "cluster_counts.csv")
    consensus_path = run_dir / "tables" / "celltypist_cluster_consensus.csv"
    consensus = pd.read_csv(consensus_path) if consensus_path.exists() else None
    registry_path = config.interpretation.evidence_registry
    registry_path = (
        registry_path if registry_path.is_absolute() else root / registry_path
    ).resolve()
    registry = load_marker_registry(registry_path)
    packet = build_evidence_packet(
        markers,
        counts,
        marker_registry=registry,
        annotation_consensus=consensus,
        top_n=config.interpretation.max_markers_per_cluster,
    )
    if config.interpretation.backend == "openai":
        model = os.getenv("OPENAI_MODEL", config.interpretation.model)
        bundle, metadata = openai_interpret(
            packet,
            model=model,
            api_key_env=config.interpretation.api_key_env,
        )
    else:
        bundle = deterministic_interpret(packet)
        metadata = {
            "backend": "deterministic",
            "validated": True,
            "external_api_called": False,
        }
    write_interpretation_artifacts(
        run_dir, packet, bundle, metadata=metadata
    )


def _run_facs_benchmark(root: Path, config: Any, run_dir: Path) -> tuple[Any, Any, Any]:
    from .benchmark import (
        build_facs_cohort,
        fetch_all_facs_datasets,
        load_label_map,
        run_facs_benchmark,
        write_benchmark_artifacts,
    )

    benchmark_config = config.benchmark
    if not benchmark_config.enabled or benchmark_config.protocol != "facs_5class":
        raise ValueError("Configuration does not enable the facs_5class benchmark")
    data_dir = benchmark_config.data_dir
    data_dir = data_dir if data_dir.is_absolute() else root / data_dir
    label_map_path = benchmark_config.label_map_path
    if label_map_path is None:
        raise ValueError("facs_5class benchmark requires label_map_path")
    label_map_path = (
        label_map_path if label_map_path.is_absolute() else root / label_map_path
    ).resolve()
    label_map = load_label_map(label_map_path)

    fetched = fetch_all_facs_datasets(data_dir)
    if benchmark_config.fail_on_same_checksum:
        digest_to_datasets: dict[str, list[str]] = {}
        for dataset_id, dataset in fetched.items():
            digest_to_datasets.setdefault(dataset.sha256, []).append(dataset_id)
        duplicates = {
            digest: dataset_ids
            for digest, dataset_ids in digest_to_datasets.items()
            if len(dataset_ids) > 1
        }
        pbmc_archive = root / "data/raw/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
        if pbmc_archive.is_file():
            from .provenance import sha256_file

            pbmc_sha = sha256_file(pbmc_archive)
            if pbmc_sha in digest_to_datasets:
                duplicates[pbmc_sha] = [
                    "pbmc3k",
                    *digest_to_datasets[pbmc_sha],
                ]
        if duplicates:
            raise ValueError(
                "benchmark source archives must be distinct from one another and "
                f"from PBMC3k: {duplicates}"
            )
    cohort = build_facs_cohort(
        {dataset_id: dataset.adata for dataset_id, dataset in fetched.items()},
        max_cells_per_class=benchmark_config.max_cells_per_class,
        random_seed=benchmark_config.random_seed,
    )
    model_dir = config.annotation.model_dir or Path("data/celltypist")
    model_dir = model_dir if model_dir.is_absolute() else root / model_dir
    result = run_facs_benchmark(
        cohort,
        model_dir,
        label_map=label_map,
        model_url=config.annotation.model_url,
        expected_model_sha256=config.annotation.model_sha256,
        target_sum=config.analysis.target_sum,
    )
    output_dir = run_dir / "benchmark"
    write_benchmark_artifacts(
        result,
        output_dir,
        cohort=cohort,
        fetched_datasets=fetched,
        label_map=label_map,
    )
    (output_dir / "benchmark_config.resolved.yaml").write_text(
        yaml.safe_dump(config.model_dump(mode="json"), sort_keys=False),
        encoding="utf-8",
    )
    return result, fetched, cohort


def run_all(
    project_root: Path,
    config_path: Path,
    *,
    run_dir_value: str | None = None,
) -> Path:
    from .annotation import annotate_cells
    from .config import load_config
    from .data import fetch_pbmc3k
    from .pipeline import run_analysis

    config = load_config(config_path)
    run_dir = _resolve_run_dir(project_root, config.output_dir, run_dir_value)
    run_dir.mkdir(parents=True, exist_ok=True)
    for name in ("logs", "tables", "plots", "data", "annotation", "provenance"):
        (run_dir / name).mkdir(parents=True, exist_ok=True)
    _write_resolved_config(run_dir, config)

    fetched = _run_stage(run_dir, "fetch", lambda: fetch_pbmc3k(project_root))
    shutil.copy2(
        fetched.manifest_path,
        run_dir / "provenance" / "pbmc3k_source_manifest.json",
    )
    adata, _, summary = _run_stage(
        run_dir,
        "analysis",
        lambda: run_analysis(
            fetched.adata,
            qc_config=config.qc,
            analysis_config=config.analysis,
            run_dir=run_dir,
        ),
    )

    if config.annotation.enabled:
        model_dir = config.annotation.model_dir or Path("data/celltypist")
        model_dir = model_dir if model_dir.is_absolute() else project_root / model_dir
        annotation_result = _run_stage(
            run_dir,
            "annotation",
            lambda: annotate_cells(
                adata,
                model_dir,
                model_url=config.annotation.model_url,
                expected_sha256=config.annotation.model_sha256,
                target_sum=config.analysis.target_sum,
            ),
        )
        _save_annotation_artifacts(adata, annotation_result, run_dir)
    else:
        _set_stage(run_dir, "annotation", "skipped", detail="disabled by configuration")

    if config.interpretation.enabled:
        _run_stage(
            run_dir,
            "interpretation",
            lambda: _run_interpretation(project_root, config, run_dir),
        )
    else:
        _set_stage(
            run_dir, "interpretation", "skipped", detail="disabled by configuration"
        )

    _record_execution(
        run_dir,
        project_root,
        operation="all",
        config_path=config_path,
    )
    _refresh_manifest(
        project_root,
        run_dir,
        config_path=config_path,
        random_seed=config.analysis.random_seed,
    )
    pointer = project_root / "results" / "LATEST_RUN.txt"
    pointer.parent.mkdir(parents=True, exist_ok=True)
    pointer.write_text(str(run_dir.relative_to(project_root)) + "\n", encoding="utf-8")
    return run_dir


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cellpilot", description="Reproducible CellPilot AI workflow"
    )
    parser.add_argument("--project-root", default=".", help="CellPilot project root")
    subparsers = parser.add_subparsers(dest="command", required=True)

    for name in ("fetch", "all"):
        command = subparsers.add_parser(name)
        command.add_argument("--config", default="configs/pbmc3k.yaml")
        if name == "all":
            command.add_argument("--run-dir")

    verify = subparsers.add_parser("verify")
    verify.add_argument("--run-dir", required=True)

    interpret = subparsers.add_parser("interpret")
    interpret.add_argument("--config", default="configs/pbmc3k.yaml")
    interpret.add_argument("--run-dir", required=True)

    benchmark = subparsers.add_parser("benchmark")
    benchmark.add_argument("--config", default="configs/benchmark.yaml")
    benchmark.add_argument("--run-dir", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    from .runtime import configure_runtime

    parser = _build_parser()
    args = parser.parse_args(argv)
    root = _resolve_project_root(args.project_root)
    configure_runtime(root)
    try:
        if args.command == "fetch":
            from .data import fetch_pbmc3k

            config_path = _resolve_config(root, args.config)
            fetched = fetch_pbmc3k(root)
            print(
                json.dumps(
                    {
                        "config": str(config_path),
                        "shape": list(fetched.adata.shape),
                        "archive": str(fetched.archive_path),
                        "sha256": fetched.sha256,
                    },
                    indent=2,
                )
            )
        elif args.command == "all":
            config_path = _resolve_config(root, args.config)
            run_dir = run_all(root, config_path, run_dir_value=args.run_dir)
            print(json.dumps({"run_dir": str(run_dir), "status": "completed"}, indent=2))
        elif args.command == "verify":
            from .provenance import verify_manifest

            run_dir = Path(args.run_dir)
            run_dir = run_dir if run_dir.is_absolute() else root / run_dir
            report = verify_manifest(run_dir)
            print(report.model_dump_json(indent=2))
            return 0 if report.valid else 1
        elif args.command == "interpret":
            from .config import load_config

            config_path = _resolve_config(root, args.config)
            config = load_config(config_path)
            run_dir = Path(args.run_dir)
            run_dir = run_dir if run_dir.is_absolute() else root / run_dir
            run_dir = run_dir.resolve()
            _run_stage(
                run_dir,
                "interpretation",
                lambda: _run_interpretation(root, config, run_dir),
            )
            _record_execution(
                run_dir,
                root,
                operation="interpret",
                config_path=config_path,
            )
            _refresh_manifest(
                root,
                run_dir,
                config_path=config_path,
                random_seed=config.analysis.random_seed,
            )
            print(json.dumps({"run_dir": str(run_dir), "status": "completed"}, indent=2))
        elif args.command == "benchmark":
            from .config import load_config

            config_path = _resolve_config(root, args.config)
            config = load_config(config_path)
            run_dir = Path(args.run_dir)
            run_dir = (run_dir if run_dir.is_absolute() else root / run_dir).resolve()
            if not run_dir.is_relative_to(root) or not run_dir.is_dir():
                raise ValueError(
                    f"Benchmark run directory must exist inside project: {run_dir}"
                )
            result, fetched, _ = _run_stage(
                run_dir,
                "benchmark",
                lambda: _run_facs_benchmark(root, config, run_dir),
            )
            _record_execution(
                run_dir,
                root,
                operation="benchmark",
                config_path=config_path,
            )
            _refresh_manifest(
                root,
                run_dir,
                config_path=config_path,
                random_seed=config.benchmark.random_seed,
            )
            print(
                json.dumps(
                    {
                        "run_dir": str(run_dir),
                        "status": "completed",
                        "metrics": result.metrics,
                    },
                    indent=2,
                )
            )
        return 0
    except Exception as exc:
        print(f"cellpilot: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
