"""Serializable schemas shared by CellPilot pipeline stages.

The models in this module deliberately contain only audit metadata. Scientific
configuration belongs in :mod:`cellpilot.config`.
"""

from __future__ import annotations

from datetime import datetime
from enum import Enum
from pathlib import Path, PurePosixPath
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator


SHA256_PATTERN = r"^[0-9a-f]{64}$"


class StrictSchema(BaseModel):
    """Base class that rejects unknown fields in persisted metadata."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class StageStatus(str, Enum):
    """A stage status that does not conflate absence with success."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    SKIPPED = "skipped"
    BLOCKED = "blocked"
    FAILED = "failed"


class FileDigest(StrictSchema):
    """Digest metadata for an input file, whose path may be external."""

    path: str
    sha256: str = Field(pattern=SHA256_PATTERN)
    size_bytes: int = Field(ge=0)


class OutputArtifact(StrictSchema):
    """Digest metadata for a file stored below a run directory."""

    path: str
    sha256: str = Field(pattern=SHA256_PATTERN)
    size_bytes: int = Field(ge=0)

    @field_validator("path")
    @classmethod
    def output_path_must_be_relative(cls, value: str) -> str:
        normalized = value.replace("\\", "/")
        path = PurePosixPath(normalized)
        if path.is_absolute() or ".." in path.parts or normalized in {"", "."}:
            raise ValueError("output artifact paths must be relative to the run directory")
        return path.as_posix()


class GitState(StrictSchema):
    """Best-effort state of the source repository at run time."""

    available: bool
    commit: str | None = None
    branch: str | None = None
    dirty: bool | None = None
    error: str | None = None


class RuntimeInfo(StrictSchema):
    """Runtime details needed to diagnose reproducibility failures."""

    python_version: str
    python_executable: str
    platform: str
    implementation: str
    packages: dict[str, str]


class RunManifest(StrictSchema):
    """Versioned manifest for a single analysis run."""

    schema_version: Literal["1.0"] = "1.0"
    run_id: str
    created_at_utc: datetime
    status: StageStatus
    command: list[str]
    random_seed: int | None = Field(default=None, ge=0)
    config_sha256: str | None = Field(default=None, pattern=SHA256_PATTERN)
    git: GitState
    runtime: RuntimeInfo
    inputs: list[FileDigest] = Field(default_factory=list)
    outputs: list[OutputArtifact] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    @field_validator("created_at_utc")
    @classmethod
    def timestamp_must_be_timezone_aware(cls, value: datetime) -> datetime:
        if value.tzinfo is None or value.utcoffset() is None:
            raise ValueError("created_at_utc must include timezone information")
        return value


class VerificationIssue(StrictSchema):
    """One mismatch found while checking a persisted run."""

    path: str
    kind: Literal[
        "missing",
        "not_a_file",
        "outside_run_directory",
        "manifest_self_reference",
        "size_mismatch",
        "sha256_mismatch",
    ]
    expected: str | int | None = None
    actual: str | int | None = None


class VerificationReport(StrictSchema):
    """Result of checking every input and output recorded in a manifest."""

    manifest_path: Path
    valid: bool
    checked_inputs: int = Field(default=0, ge=0)
    checked_outputs: int = Field(ge=0)
    issues: list[VerificationIssue] = Field(default_factory=list)
