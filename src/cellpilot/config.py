"""Strict YAML configuration models for CellPilot AI."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal
from urllib.parse import urlparse

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class StrictConfig(BaseModel):
    """Configuration base class that fails on misspelled/unknown keys."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class DatasetConfig(StrictConfig):
    """Source and integrity metadata for one expression dataset."""

    name: str = Field(default="pbmc3k", min_length=1)
    source: Literal["scanpy", "url", "local"] = "scanpy"
    loader: str | None = "pbmc3k"
    url: str | None = None
    path: Path | None = None
    sha256: str | None = None
    format: Literal["h5ad", "10x_mtx", "mtx", "csv", "tsv"] = "h5ad"
    description: str | None = None
    license: str | None = None

    @field_validator("url")
    @classmethod
    def validate_http_url(cls, value: str | None) -> str | None:
        if value is None:
            return None
        parsed = urlparse(value)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            raise ValueError("url must be an absolute http(s) URL")
        return value

    @field_validator("sha256")
    @classmethod
    def normalize_sha256(cls, value: str | None) -> str | None:
        if value is None:
            return None
        normalized = value.strip().lower()
        if len(normalized) != 64 or any(c not in "0123456789abcdef" for c in normalized):
            raise ValueError("sha256 must contain exactly 64 hexadecimal characters")
        return normalized

    @model_validator(mode="after")
    def validate_source_fields(self) -> "DatasetConfig":
        if self.source == "scanpy" and not self.loader:
            raise ValueError("loader is required when source='scanpy'")
        if self.source == "url" and not self.url:
            raise ValueError("url is required when source='url'")
        if self.source == "local" and self.path is None:
            raise ValueError("path is required when source='local'")
        return self


class QCConfig(StrictConfig):
    """Cell- and gene-filtering thresholds applied to raw counts."""

    min_genes: int = Field(default=200, ge=0)
    max_genes: int | None = Field(default=2_500, gt=0)
    min_counts: int | None = Field(default=None, ge=0)
    max_counts: int | None = Field(default=None, gt=0)
    max_pct_mt: float | None = Field(default=5.0, ge=0, le=100)
    min_cells_per_gene: int = Field(default=3, ge=1)
    mito_prefix: tuple[str, ...] = ("MT-",)

    @field_validator("mito_prefix")
    @classmethod
    def prefixes_must_be_nonempty(cls, value: tuple[str, ...]) -> tuple[str, ...]:
        if not value or any(not prefix for prefix in value):
            raise ValueError("mito_prefix must contain at least one non-empty prefix")
        return value

    @model_validator(mode="after")
    def validate_threshold_order(self) -> "QCConfig":
        if self.max_genes is not None and self.max_genes < self.min_genes:
            raise ValueError("max_genes must be greater than or equal to min_genes")
        if (
            self.min_counts is not None
            and self.max_counts is not None
            and self.max_counts < self.min_counts
        ):
            raise ValueError("max_counts must be greater than or equal to min_counts")
        return self


class AnalysisConfig(StrictConfig):
    """Deterministic normalization, embedding, clustering, and marker settings."""

    target_sum: float = Field(default=10_000.0, gt=0)
    n_top_genes: int = Field(default=2_000, gt=0)
    hvg_flavor: Literal["seurat", "seurat_v3", "cell_ranger"] = "seurat"
    n_pcs: int = Field(default=40, gt=1)
    n_neighbors: int = Field(default=10, gt=1)
    umap_min_dist: float = Field(default=0.5, ge=0, le=1)
    leiden_resolution: float = Field(default=0.5, gt=0)
    marker_method: Literal["wilcoxon", "t-test", "t-test_overestim_var", "logreg"] = (
        "wilcoxon"
    )
    marker_n_genes: int = Field(default=100, gt=0)
    random_seed: int = Field(default=42, ge=0)

    @property
    def marker_top_n(self) -> int:
        """Compatibility alias for pipeline code; YAML uses ``marker_n_genes``."""

        return self.marker_n_genes


class AnnotationConfig(StrictConfig):
    """CellTypist adapter settings.

    Model download metadata are optional so an unavailable remote model can be
    represented honestly as blocked instead of preventing deterministic stages.
    """

    enabled: bool = True
    method: Literal["celltypist"] = "celltypist"
    model: Literal["Healthy_COVID19_PBMC.pkl"] = "Healthy_COVID19_PBMC.pkl"
    model_dir: Path | None = None
    model_url: str | None = (
        "https://celltypist.cog.sanger.ac.uk/models/COVID19_PBMC_Haniffa/v1/"
        "Healthy_COVID19_PBMC.pkl"
    )
    model_sha256: str | None = (
        "1eee68d23b4ac0f6999a0a92d6c950d671f31fcdb003ce0760c580b0bd6a3b79"
    )
    majority_voting: Literal[False] = False
    strict: Literal[True] = True
    min_confidence: Literal[0.0] = 0.0

    @field_validator("model_url")
    @classmethod
    def validate_model_url(cls, value: str | None) -> str | None:
        if value is None:
            return None
        parsed = urlparse(value)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            raise ValueError("model_url must be an absolute http(s) URL")
        expected = (
            "https://celltypist.cog.sanger.ac.uk/models/COVID19_PBMC_Haniffa/v1/"
            "Healthy_COVID19_PBMC.pkl"
        )
        if value != expected:
            raise ValueError("model_url must identify the frozen CellTypist artifact")
        return value

    @field_validator("model_sha256")
    @classmethod
    def normalize_model_sha256(cls, value: str | None) -> str | None:
        normalized = DatasetConfig.normalize_sha256(value)
        expected = "1eee68d23b4ac0f6999a0a92d6c950d671f31fcdb003ce0760c580b0bd6a3b79"
        if normalized is not None and normalized != expected:
            raise ValueError("model_sha256 must match the frozen CellTypist artifact")
        return normalized

    @property
    def model_name(self) -> str:
        """Compatibility alias used by the CellTypist adapter."""

        return self.model


class InterpretationConfig(StrictConfig):
    """Evidence-grounded deterministic or optional LLM interpretation."""

    enabled: bool = True
    backend: Literal["deterministic", "openai"] = "deterministic"
    model: str = Field(default="gpt-5-mini", min_length=1)
    api_key_env: str = Field(default="OPENAI_API_KEY", min_length=1)
    evidence_registry: Path = Path("src/cellpilot/resources/marker_evidence.yaml")
    max_markers_per_cluster: int = Field(default=20, gt=0)
    require_citations: Literal[True] = True
    reject_unsupported: Literal[True] = True


class BenchmarkConfig(StrictConfig):
    """Independent labeled-dataset benchmark settings."""

    enabled: bool = False
    protocol: Literal["single_dataset", "facs_5class"] = "single_dataset"
    dataset: DatasetConfig | None = None
    label_column: str | None = None
    label_map_path: Path | None = None
    data_dir: Path = Path("data/benchmark")
    max_cells_per_class: int | None = Field(default=None, gt=0)
    bootstrap_iterations: int = Field(default=1_000, ge=0)
    confidence_level: float = Field(default=0.95, gt=0, lt=1)
    random_seed: int = Field(default=42, ge=0)
    fail_on_same_checksum: bool = True

    @model_validator(mode="after")
    def require_enabled_fields(self) -> "BenchmarkConfig":
        if self.enabled and self.protocol == "single_dataset" and self.dataset is None:
            raise ValueError("dataset is required when benchmark.enabled=true")
        if self.enabled and self.protocol == "single_dataset" and not self.label_column:
            raise ValueError("label_column is required when benchmark.enabled=true")
        if self.enabled and self.protocol == "facs_5class" and self.label_map_path is None:
            raise ValueError("label_map_path is required for the facs_5class benchmark")
        if self.protocol == "facs_5class" and self.bootstrap_iterations != 0:
            raise ValueError(
                "bootstrap_iterations must be 0 for facs_5class because cell-level "
                "resampling cannot supply missing biological replication"
            )
        return self


class ProjectConfig(StrictConfig):
    """Top-level CellPilot project configuration."""

    project_name: str = Field(default="cellpilot-ai", min_length=1)
    output_dir: Path = Path("results")
    dataset: DatasetConfig = Field(default_factory=DatasetConfig)
    qc: QCConfig = Field(default_factory=QCConfig)
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)
    interpretation: InterpretationConfig = Field(default_factory=InterpretationConfig)
    benchmark: BenchmarkConfig = Field(default_factory=BenchmarkConfig)


# A descriptive alias for callers that prefer the application name.
CellPilotConfig = ProjectConfig


def load_config(path: str | Path) -> ProjectConfig:
    """Load and strictly validate one UTF-8 YAML configuration file.

    Relative paths are intentionally preserved. CLI callers should resolve them
    against their explicitly selected project root, not the process home or an
    implicit global directory.
    """

    config_path = Path(path)
    try:
        payload: Any = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ValueError(f"invalid YAML in {config_path}: {exc}") from exc

    if payload is None:
        payload = {}
    if not isinstance(payload, dict):
        raise ValueError(f"configuration root in {config_path} must be a mapping")
    return ProjectConfig.model_validate(payload)
