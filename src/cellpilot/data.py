"""Auditable public-data acquisition and raw-count validation."""

from __future__ import annotations

import hashlib
import json
import shutil
import tarfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np
from scipy import sparse

PBMC3K_URL = (
    "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
    "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)
PBMC3K_SHA256 = "847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5"


@dataclass(frozen=True)
class FetchedDataset:
    adata: Any
    archive_path: Path
    matrix_dir: Path
    manifest_path: Path
    sha256: str


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(
    url: str,
    destination: Path,
    *,
    expected_sha256: str | None = None,
    timeout_seconds: int = 120,
) -> str:
    """Stream a URL to a workspace path and verify its digest before promotion."""

    import requests

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    digest = hashlib.sha256()
    try:
        with requests.get(url, stream=True, timeout=timeout_seconds) as response:
            response.raise_for_status()
            with partial.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        digest.update(chunk)
                        handle.write(chunk)
        observed = digest.hexdigest()
        if expected_sha256 and observed.lower() != expected_sha256.lower():
            raise ValueError(
                f"Checksum mismatch for {url}: expected {expected_sha256}, observed {observed}"
            )
        partial.replace(destination)
        return observed
    finally:
        if partial.exists():
            partial.unlink()


def _archive_target(member: tarfile.TarInfo, destination: Path) -> Path:
    normalized = member.name.replace("\\", "/")
    relative = PurePosixPath(normalized)
    if (
        normalized in {"", "."}
        or relative.is_absolute()
        or ".." in relative.parts
        or not relative.parts
    ):
        raise ValueError(f"Unsafe path in archive: {member.name}")
    target = (destination / Path(*relative.parts)).resolve()
    resolved_destination = destination.resolve()
    if target != resolved_destination and resolved_destination not in target.parents:
        raise ValueError(f"Unsafe path in archive: {member.name}")
    return target


def _safe_extract_tar(archive: Path, destination: Path) -> None:
    """Extract only regular files/directories after rejecting unsafe members."""

    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle.getmembers():
            _archive_target(member, destination)
            if not (member.isdir() or member.isfile()):
                raise ValueError(
                    f"Archive member is not a regular file/directory: {member.name}"
                )
        bundle.extractall(destination, filter="data")


def verify_extracted_tar(archive: Path, destination: Path) -> dict[str, str]:
    """Compare every extracted regular file byte-for-byte with its tar member."""

    observed: dict[str, str] = {}
    with tarfile.open(archive, "r:*") as bundle:
        for member in bundle.getmembers():
            target = _archive_target(member, destination)
            if member.isdir():
                if not target.is_dir():
                    raise ValueError(f"Extracted directory is missing: {member.name}")
                continue
            if not member.isfile():
                raise ValueError(
                    f"Archive member is not a regular file/directory: {member.name}"
                )
            if not target.is_file() or target.is_symlink():
                raise ValueError(f"Extracted file is missing or unsafe: {member.name}")
            source = bundle.extractfile(member)
            if source is None:
                raise ValueError(f"Could not read archive member: {member.name}")
            archive_digest = hashlib.sha256()
            with source:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    archive_digest.update(chunk)
            extracted_digest = sha256_file(target)
            expected = archive_digest.hexdigest()
            if target.stat().st_size != member.size or extracted_digest != expected:
                raise ValueError(
                    f"Extracted file does not match verified archive: {member.name}"
                )
            observed[member.name] = extracted_digest
    return observed


def _find_10x_matrix_dir(root: Path) -> Path:
    candidates = sorted(
        {
            path.parent
            for pattern in ("matrix.mtx", "matrix.mtx.gz")
            for path in root.rglob(pattern)
        }
    )
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected exactly one 10x matrix directory below {root}, found {len(candidates)}"
        )
    return candidates[0]


def _matrix_values(matrix: Any) -> np.ndarray:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    return np.asarray(values)


def validate_raw_counts(adata: Any, *, atol: float = 1e-6) -> dict[str, float | int]:
    """Assert that an AnnData matrix is finite, non-negative, integer-like counts."""

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Raw count matrix is empty")
    values = _matrix_values(adata.X)
    if not np.isfinite(values).all():
        raise ValueError("Raw count matrix contains non-finite values")
    if values.size and float(values.min()) < 0:
        raise ValueError("Raw count matrix contains negative values")
    if values.size and not np.allclose(values, np.rint(values), atol=atol, rtol=0):
        raise ValueError("Raw count matrix is not integer-like")
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "nonzero_values": int(values.size),
        "max_count": float(values.max()) if values.size else 0.0,
    }


def fetch_pbmc3k(project_root: Path, *, force: bool = False) -> FetchedDataset:
    """Download the official 10x PBMC 3k filtered matrix and load raw counts."""

    import scanpy as sc

    raw_dir = project_root.resolve() / "data" / "raw" / "pbmc3k"
    archive = raw_dir / "pbmc3k_filtered_gene_bc_matrices.tar.gz"
    extracted = raw_dir / "extracted"
    manifest_path = raw_dir / "source_manifest.json"

    if force and extracted.exists():
        shutil.rmtree(extracted)
    if force or not archive.exists():
        observed_sha = download_file(
            PBMC3K_URL, archive, expected_sha256=PBMC3K_SHA256
        )
    else:
        observed_sha = sha256_file(archive)
        if observed_sha != PBMC3K_SHA256:
            raise ValueError(
                f"Cached PBMC3k archive checksum mismatch: {observed_sha}"
            )
    if not extracted.exists():
        _safe_extract_tar(archive, extracted)
    extracted_digests = verify_extracted_tar(archive, extracted)

    matrix_dir = _find_10x_matrix_dir(extracted)
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)
    adata.var["gene_symbol_original"] = adata.var_names.astype(str)
    adata.var_names_make_unique()
    validation = validate_raw_counts(adata)
    manifest = {
        "dataset": "10x PBMC 3k",
        "source_url": PBMC3K_URL,
        "archive_sha256_expected": PBMC3K_SHA256,
        "archive_sha256_observed": observed_sha,
        "license": "CC BY 4.0 (as stated by the 10x Genomics dataset page)",
        "input_kind": "Cell Ranger filtered gene-barcode UMI count matrix; not FASTQ or unfiltered droplets",
        "matrix_dir": str(matrix_dir.relative_to(project_root.resolve())),
        "extracted_member_sha256": extracted_digests,
        **validation,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return FetchedDataset(
        adata=adata,
        archive_path=archive,
        matrix_dir=matrix_dir,
        manifest_path=manifest_path,
        sha256=observed_sha,
    )
