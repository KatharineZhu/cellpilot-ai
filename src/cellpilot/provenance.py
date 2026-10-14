"""Reproducibility metadata and artifact-integrity verification."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Iterable, Sequence

from .schemas import (
    FileDigest,
    GitState,
    OutputArtifact,
    RunManifest,
    RuntimeInfo,
    StageStatus,
    VerificationIssue,
    VerificationReport,
)


DEFAULT_PACKAGES = (
    "anndata",
    "certifi",
    "celltypist",
    "h5py",
    "igraph",
    "leidenalg",
    "matplotlib",
    "numba",
    "numpy",
    "openai",
    "pandas",
    "pydantic",
    "scanpy",
    "scikit-learn",
    "scipy",
    "streamlit",
    "umap-learn",
)
MANIFEST_FILENAME = "manifest.json"


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of a file without loading it into memory."""

    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    file_path = Path(path)
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run_git(project_root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", "-C", str(project_root), *args],
        check=False,
        capture_output=True,
        text=True,
        timeout=5,
    )


def collect_git_state(project_root: str | Path) -> GitState:
    """Collect repository state without modifying the repository."""

    root = Path(project_root).resolve()
    try:
        inside = _run_git(root, "rev-parse", "--is-inside-work-tree")
    except (FileNotFoundError, subprocess.SubprocessError) as exc:
        return GitState(available=False, error=str(exc))
    if inside.returncode != 0 or inside.stdout.strip() != "true":
        error = inside.stderr.strip() or "not a Git worktree"
        return GitState(available=False, error=error)

    commit_result = _run_git(root, "rev-parse", "HEAD")
    branch_result = _run_git(root, "symbolic-ref", "--quiet", "--short", "HEAD")
    status_result = _run_git(root, "status", "--porcelain", "--untracked-files=normal")
    commit = commit_result.stdout.strip() if commit_result.returncode == 0 else None
    branch = branch_result.stdout.strip() if branch_result.returncode == 0 else None
    dirty = bool(status_result.stdout.strip()) if status_result.returncode == 0 else None
    errors = [
        result.stderr.strip()
        for result in (commit_result, branch_result, status_result)
        if result.returncode != 0 and result.stderr.strip()
    ]
    return GitState(
        available=True,
        commit=commit,
        branch=branch,
        dirty=dirty,
        error="; ".join(errors) or None,
    )


def collect_package_versions(packages: Iterable[str] = DEFAULT_PACKAGES) -> dict[str, str]:
    """Return installed versions, explicitly marking absent distributions."""

    versions: dict[str, str] = {}
    for package in sorted(set(packages)):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def collect_runtime(packages: Iterable[str] = DEFAULT_PACKAGES) -> RuntimeInfo:
    """Capture Python, operating system, and relevant dependency versions."""

    return RuntimeInfo(
        python_version=platform.python_version(),
        python_executable=sys.executable,
        platform=platform.platform(),
        implementation=platform.python_implementation(),
        packages=collect_package_versions(packages),
    )


def _contained_file(path: Path, root: Path) -> Path:
    resolved_root = root.resolve()
    resolved_path = path.resolve()
    if not resolved_path.is_relative_to(resolved_root):
        raise ValueError(f"output path is outside run directory: {path}")
    return resolved_path


def _expand_output_paths(run_dir: Path, outputs: Sequence[str | Path] | None) -> list[Path]:
    if outputs is None:
        candidates = list(run_dir.rglob("*")) if run_dir.exists() else []
    else:
        candidates = []
        for value in outputs:
            candidate = Path(value)
            if not candidate.is_absolute():
                candidate = run_dir / candidate
            candidate = _contained_file(candidate, run_dir)
            if candidate.is_dir():
                candidates.extend(candidate.rglob("*"))
            else:
                candidates.append(candidate)

    manifest_path = (run_dir / MANIFEST_FILENAME).resolve()
    temporary_manifest_path = (run_dir / f".{MANIFEST_FILENAME}.tmp").resolve()
    files: set[Path] = set()
    for candidate in candidates:
        resolved = _contained_file(candidate, run_dir)
        if resolved in {manifest_path, temporary_manifest_path}:
            continue
        if not resolved.exists():
            raise FileNotFoundError(resolved)
        if not resolved.is_file():
            continue
        files.add(resolved)
    return sorted(files, key=lambda item: item.relative_to(run_dir.resolve()).as_posix())


def _digest_input(path: str | Path) -> FileDigest:
    file_path = Path(path).resolve()
    if not file_path.is_file():
        raise FileNotFoundError(file_path)
    return FileDigest(
        path=str(file_path),
        sha256=sha256_file(file_path),
        size_bytes=file_path.stat().st_size,
    )


def build_manifest(
    run_dir: str | Path,
    *,
    project_root: str | Path,
    outputs: Sequence[str | Path] | None = None,
    inputs: Sequence[str | Path] = (),
    config_path: str | Path | None = None,
    command: Sequence[str] | None = None,
    packages: Iterable[str] = DEFAULT_PACKAGES,
    random_seed: int | None = None,
    status: StageStatus | str = StageStatus.COMPLETED,
    warnings: Sequence[str] = (),
    created_at_utc: datetime | None = None,
) -> RunManifest:
    """Construct a manifest from files already present in ``run_dir``.

    The manifest itself is always excluded, including when ``outputs`` is left
    as ``None`` and files are discovered recursively.
    """

    resolved_run_dir = Path(run_dir).resolve()
    output_records = [
        OutputArtifact(
            path=path.relative_to(resolved_run_dir).as_posix(),
            sha256=sha256_file(path),
            size_bytes=path.stat().st_size,
        )
        for path in _expand_output_paths(resolved_run_dir, outputs)
    ]
    input_records = [_digest_input(path) for path in inputs]
    config_sha256 = None
    if config_path is not None:
        config_file = Path(config_path)
        if not config_file.is_file():
            raise FileNotFoundError(config_file)
        config_sha256 = sha256_file(config_file)

    timestamp = created_at_utc or datetime.now(timezone.utc)
    return RunManifest(
        run_id=resolved_run_dir.name,
        created_at_utc=timestamp,
        status=status,
        command=list(command) if command is not None else list(sys.argv),
        random_seed=random_seed,
        config_sha256=config_sha256,
        git=collect_git_state(project_root),
        runtime=collect_runtime(packages),
        inputs=input_records,
        outputs=output_records,
        warnings=list(warnings),
    )


def write_manifest(
    run_dir: str | Path,
    *,
    project_root: str | Path,
    outputs: Sequence[str | Path] | None = None,
    inputs: Sequence[str | Path] = (),
    config_path: str | Path | None = None,
    command: Sequence[str] | None = None,
    packages: Iterable[str] = DEFAULT_PACKAGES,
    random_seed: int | None = None,
    status: StageStatus | str = StageStatus.COMPLETED,
    warnings: Sequence[str] = (),
    created_at_utc: datetime | None = None,
) -> Path:
    """Write ``manifest.json`` atomically inside the explicitly passed run dir."""

    resolved_run_dir = Path(run_dir).resolve()
    resolved_run_dir.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest(
        resolved_run_dir,
        project_root=project_root,
        outputs=outputs,
        inputs=inputs,
        config_path=config_path,
        command=command,
        packages=packages,
        random_seed=random_seed,
        status=status,
        warnings=warnings,
        created_at_utc=created_at_utc,
    )
    manifest_path = resolved_run_dir / MANIFEST_FILENAME
    temporary_path = resolved_run_dir / f".{MANIFEST_FILENAME}.tmp"
    serialized = manifest.model_dump_json(indent=2) + "\n"
    temporary_path.write_text(serialized, encoding="utf-8")
    temporary_path.replace(manifest_path)
    return manifest_path


def load_manifest(path: str | Path) -> RunManifest:
    """Load a manifest from a file or a run directory."""

    manifest_path = Path(path)
    if manifest_path.is_dir():
        manifest_path = manifest_path / MANIFEST_FILENAME
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    return RunManifest.model_validate(payload)


def verify_manifest(path: str | Path) -> VerificationReport:
    """Verify recorded input/output sizes and hashes without hashing the manifest."""

    supplied_path = Path(path)
    manifest_path = supplied_path / MANIFEST_FILENAME if supplied_path.is_dir() else supplied_path
    manifest_path = manifest_path.resolve()
    run_dir = manifest_path.parent
    manifest = load_manifest(manifest_path)
    issues: list[VerificationIssue] = []
    checked_inputs = 0
    checked_outputs = 0

    for artifact in manifest.inputs:
        input_path = Path(artifact.path).expanduser()
        if not input_path.exists():
            issues.append(VerificationIssue(path=artifact.path, kind="missing"))
            continue
        if not input_path.is_file():
            issues.append(VerificationIssue(path=artifact.path, kind="not_a_file"))
            continue

        checked_inputs += 1
        actual_size = input_path.stat().st_size
        if actual_size != artifact.size_bytes:
            issues.append(
                VerificationIssue(
                    path=artifact.path,
                    kind="size_mismatch",
                    expected=artifact.size_bytes,
                    actual=actual_size,
                )
            )
        actual_sha256 = sha256_file(input_path)
        if actual_sha256 != artifact.sha256:
            issues.append(
                VerificationIssue(
                    path=artifact.path,
                    kind="sha256_mismatch",
                    expected=artifact.sha256,
                    actual=actual_sha256,
                )
            )

    for artifact in manifest.outputs:
        artifact_path = run_dir / artifact.path
        try:
            resolved = _contained_file(artifact_path, run_dir)
        except ValueError:
            issues.append(
                VerificationIssue(path=artifact.path, kind="outside_run_directory")
            )
            continue
        if resolved == manifest_path:
            issues.append(VerificationIssue(path=artifact.path, kind="manifest_self_reference"))
            continue
        if not resolved.exists():
            issues.append(VerificationIssue(path=artifact.path, kind="missing"))
            continue
        if not resolved.is_file():
            issues.append(VerificationIssue(path=artifact.path, kind="not_a_file"))
            continue

        checked_outputs += 1
        actual_size = resolved.stat().st_size
        if actual_size != artifact.size_bytes:
            issues.append(
                VerificationIssue(
                    path=artifact.path,
                    kind="size_mismatch",
                    expected=artifact.size_bytes,
                    actual=actual_size,
                )
            )
        actual_sha256 = sha256_file(resolved)
        if actual_sha256 != artifact.sha256:
            issues.append(
                VerificationIssue(
                    path=artifact.path,
                    kind="sha256_mismatch",
                    expected=artifact.sha256,
                    actual=actual_sha256,
                )
            )

    return VerificationReport(
        manifest_path=manifest_path,
        valid=not issues,
        checked_inputs=checked_inputs,
        checked_outputs=checked_outputs,
        issues=issues,
    )
