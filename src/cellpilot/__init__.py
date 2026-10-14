"""CellPilot AI: reproducible single-cell analysis with auditable outputs."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from .config import CellPilotConfig, ProjectConfig, load_config

try:
    __version__ = version("cellpilot-ai")
except PackageNotFoundError:  # Source checkout before editable installation.
    __version__ = "0.1.0"

__all__ = ["CellPilotConfig", "ProjectConfig", "__version__", "load_config"]
