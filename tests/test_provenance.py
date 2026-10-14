from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from cellpilot.provenance import (
    collect_git_state,
    collect_runtime,
    sha256_file,
    verify_manifest,
    write_manifest,
)


def test_sha256_file_matches_known_digest(tmp_path: Path) -> None:
    path = tmp_path / "payload.bin"
    path.write_bytes(b"abc")

    assert sha256_file(path) == hashlib.sha256(b"abc").hexdigest()


def test_manifest_records_outputs_and_excludes_itself(tmp_path: Path) -> None:
    project_root = tmp_path / "project"
    run_dir = project_root / "results" / "run-001"
    table = run_dir / "tables" / "markers.csv"
    table.parent.mkdir(parents=True)
    table.write_text("cluster,gene\n0,MS4A1\n", encoding="utf-8")
    config = project_root / "config.yaml"
    config.write_text("project_name: test\n", encoding="utf-8")
    input_file = project_root / "input.h5ad"
    input_file.write_bytes(b"fixture")

    manifest_path = write_manifest(
        run_dir,
        project_root=project_root,
        inputs=[input_file],
        config_path=config,
        command=["cellpilot", "analyze"],
        packages=("pydantic", "package-that-does-not-exist-cellpilot"),
        random_seed=42,
        created_at_utc=datetime(2026, 1, 2, tzinfo=timezone.utc),
    )

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert [item["path"] for item in payload["outputs"]] == ["tables/markers.csv"]
    assert payload["random_seed"] == 42
    assert payload["config_sha256"] == sha256_file(config)
    assert payload["runtime"]["packages"]["package-that-does-not-exist-cellpilot"] == (
        "not-installed"
    )
    report = verify_manifest(run_dir)
    assert report.valid
    assert report.checked_inputs == 1
    assert report.checked_outputs == 1
    assert report.issues == []


def test_verify_manifest_detects_tampering(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    artifact = run_dir / "result.txt"
    artifact.parent.mkdir()
    artifact.write_text("original", encoding="utf-8")
    write_manifest(
        run_dir,
        project_root=tmp_path,
        outputs=[artifact],
        packages=(),
        command=[],
    )

    artifact.write_text("tampered and longer", encoding="utf-8")
    report = verify_manifest(run_dir)

    assert not report.valid
    assert {issue.kind for issue in report.issues} == {"size_mismatch", "sha256_mismatch"}


def test_verify_manifest_detects_input_tampering(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    source = tmp_path / "source.bin"
    source.write_bytes(b"original input")
    write_manifest(
        run_dir,
        project_root=tmp_path,
        inputs=[source],
        packages=(),
        command=[],
    )

    source.write_bytes(b"tampered input with a different size")
    report = verify_manifest(run_dir)

    assert not report.valid
    assert report.checked_inputs == 1
    assert {issue.kind for issue in report.issues} == {
        "size_mismatch",
        "sha256_mismatch",
    }


def test_write_manifest_rejects_output_outside_run_directory(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    outside = tmp_path / "outside.txt"
    outside.write_text("do not include", encoding="utf-8")

    with pytest.raises(ValueError, match="outside run directory"):
        write_manifest(
            run_dir,
            project_root=tmp_path,
            outputs=[outside],
            packages=(),
        )


def test_collect_git_state_handles_non_repository(tmp_path: Path) -> None:
    state = collect_git_state(tmp_path)

    assert not state.available
    assert state.commit is None
    assert state.error


def test_collect_runtime_is_network_free_and_marks_missing_package() -> None:
    runtime = collect_runtime(("pydantic", "package-that-does-not-exist-cellpilot"))

    assert runtime.python_version
    assert runtime.platform
    assert runtime.packages["pydantic"] != "not-installed"
    assert runtime.packages["package-that-does-not-exist-cellpilot"] == "not-installed"
