from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from cellpilot.config import (
    AnnotationConfig,
    BenchmarkConfig,
    DatasetConfig,
    ProjectConfig,
    QCConfig,
    load_config,
)


def test_load_minimal_config_uses_documented_defaults(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("project_name: test-project\n", encoding="utf-8")

    config = load_config(config_path)

    assert config.project_name == "test-project"
    assert config.dataset.name == "pbmc3k"
    assert config.dataset.source == "scanpy"
    assert config.qc.max_pct_mt == 5.0
    assert config.analysis.random_seed == 42
    assert config.output_dir == Path("results")


def test_load_full_nested_config(tmp_path: Path) -> None:
    checksum = "AB" * 32
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        f"""
dataset:
  name: remote-test
  source: url
  loader: null
  url: https://example.org/cells.h5ad
  sha256: {checksum}
  format: h5ad
qc:
  min_genes: 10
  max_genes: 1000
  max_pct_mt: 10
  min_cells_per_gene: 2
  mito_prefix: [MT-, mt-]
analysis:
  random_seed: 7
annotation:
  enabled: false
interpretation:
  backend: deterministic
benchmark:
  enabled: false
""".lstrip(),
        encoding="utf-8",
    )

    config = load_config(config_path)

    assert config.dataset.sha256 == checksum.lower()
    assert config.qc.mito_prefix == ("MT-", "mt-")
    assert config.analysis.random_seed == 7
    assert not config.annotation.enabled


def test_unknown_key_is_rejected(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("qc:\n  max_pct_mito: 5\n", encoding="utf-8")

    with pytest.raises(ValidationError, match="max_pct_mito"):
        load_config(config_path)


@pytest.mark.parametrize(
    "payload",
    [
        {"min_genes": 500, "max_genes": 100},
        {"min_counts": 1000, "max_counts": 100},
        {"max_pct_mt": 101},
        {"mito_prefix": []},
    ],
)
def test_invalid_qc_thresholds_are_rejected(payload: dict[str, object]) -> None:
    with pytest.raises(ValidationError):
        QCConfig.model_validate(payload)


def test_dataset_source_requires_matching_locator() -> None:
    with pytest.raises(ValidationError, match="url is required"):
        DatasetConfig(source="url", loader=None)
    with pytest.raises(ValidationError, match="path is required"):
        DatasetConfig(source="local", loader=None)
    with pytest.raises(ValidationError, match="loader is required"):
        DatasetConfig(source="scanpy", loader=None)


def test_enabled_benchmark_requires_dataset_and_label() -> None:
    with pytest.raises(ValidationError, match="dataset is required"):
        BenchmarkConfig(enabled=True, label_column="cell_type")

    dataset = DatasetConfig(name="reference", source="scanpy", loader="pbmc68k_reduced")
    with pytest.raises(ValidationError, match="label_column is required"):
        BenchmarkConfig(enabled=True, dataset=dataset)


def test_facs_benchmark_uses_frozen_label_map_without_cell_label_column() -> None:
    config = BenchmarkConfig(
        enabled=True,
        protocol="facs_5class",
        label_map_path=Path("src/cellpilot/resources/benchmark_label_map.yaml"),
        max_cells_per_class=500,
        bootstrap_iterations=0,
    )
    assert config.dataset is None
    assert config.label_column is None


@pytest.mark.parametrize(
    "override",
    [
        {"majority_voting": True},
        {"strict": False},
        {"min_confidence": 0.5},
        {"model": "Immune_All_Low.pkl"},
    ],
)
def test_annotation_config_rejects_unsupported_knobs(
    override: dict[str, object],
) -> None:
    with pytest.raises(ValidationError):
        AnnotationConfig.model_validate(override)


def test_facs_benchmark_rejects_cell_level_bootstrap_claim() -> None:
    with pytest.raises(ValidationError, match="bootstrap_iterations must be 0"):
        BenchmarkConfig(
            enabled=True,
            protocol="facs_5class",
            label_map_path=Path("mapping.yaml"),
            bootstrap_iterations=100,
        )


def test_project_config_rejects_unknown_top_level_field() -> None:
    with pytest.raises(ValidationError):
        ProjectConfig.model_validate({"seed": 42})


def test_non_mapping_yaml_is_rejected(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("- not\n- a\n- mapping\n", encoding="utf-8")

    with pytest.raises(ValueError, match="must be a mapping"):
        load_config(config_path)
