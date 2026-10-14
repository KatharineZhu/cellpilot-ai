from __future__ import annotations

import io
import tarfile

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from cellpilot.data import validate_raw_counts, verify_extracted_tar


def test_validate_raw_counts_accepts_sparse_integer_matrix() -> None:
    adata = ad.AnnData(sparse.csr_matrix([[0, 1, 2], [3, 0, 1]], dtype=float))
    result = validate_raw_counts(adata)
    assert result["n_cells"] == 2
    assert result["n_genes"] == 3
    assert result["max_count"] == 3


@pytest.mark.parametrize(
    "matrix, message",
    [
        (np.array([[0.0, -1.0]]), "negative"),
        (np.array([[0.0, 1.25]]), "integer-like"),
        (np.array([[0.0, np.nan]]), "non-finite"),
    ],
)
def test_validate_raw_counts_rejects_invalid_values(matrix: np.ndarray, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        validate_raw_counts(ad.AnnData(matrix))


def test_verify_extracted_tar_detects_modified_cache(tmp_path) -> None:
    archive = tmp_path / "source.tar.gz"
    payload = b"verified matrix bytes"
    with tarfile.open(archive, "w:gz") as bundle:
        member = tarfile.TarInfo("matrix/matrix.mtx")
        member.size = len(payload)
        bundle.addfile(member, io.BytesIO(payload))

    extracted = tmp_path / "extracted"
    target = extracted / "matrix" / "matrix.mtx"
    target.parent.mkdir(parents=True)
    target.write_bytes(payload)
    assert verify_extracted_tar(archive, extracted)["matrix/matrix.mtx"]

    target.write_bytes(b"tampered matrix bytes")
    with pytest.raises(ValueError, match="does not match verified archive"):
        verify_extracted_tar(archive, extracted)
