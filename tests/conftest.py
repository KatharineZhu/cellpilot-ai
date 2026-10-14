from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / "work" / "test-cache" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(PROJECT_ROOT / "work" / "test-cache" / "numba"))
os.environ.setdefault("CELLTYPIST_FOLDER", str(PROJECT_ROOT / "work" / "test-cache" / "celltypist"))

for variable in ("MPLCONFIGDIR", "NUMBA_CACHE_DIR", "CELLTYPIST_FOLDER"):
    Path(os.environ[variable]).mkdir(parents=True, exist_ok=True)
