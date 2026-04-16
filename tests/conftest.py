"""Prepend ``src/`` so local ``pytest`` resolves the in-tree ``gwaspeek`` package first."""

from __future__ import annotations

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
_src_s = str(_SRC)
if _src_s not in sys.path:
    sys.path.insert(0, _src_s)
