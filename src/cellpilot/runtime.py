"""Runtime isolation helpers.

All caches created by CellPilot are kept below the project root so a run does not
silently write scientific state into a user's home directory.
"""

from __future__ import annotations

import os
from pathlib import Path


def configure_runtime(project_root: Path) -> dict[str, Path]:
    """Configure writable, project-local cache directories.

    This function must be called before importing matplotlib, Scanpy, Numba, or
    CellTypist in command-line entry points.
    """

    root = project_root.resolve()
    paths = {
        "matplotlib": root / "work" / "cache" / "matplotlib",
        "numba": root / "work" / "cache" / "numba",
        "celltypist": root / "data" / "celltypist",
        "scanpy": root / "data" / "cache" / "scanpy",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(paths["matplotlib"]))
    os.environ.setdefault("NUMBA_CACHE_DIR", str(paths["numba"]))
    os.environ.setdefault("CELLTYPIST_FOLDER", str(paths["celltypist"]))
    return paths
