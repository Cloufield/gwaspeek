from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import pandas as pd


CHR_ALIASES: Dict[str, int] = {"X": 23, "Y": 24, "MT": 25, "M": 25}


def _normalize_chr_token(value: object) -> Optional[int]:
    if pd.isna(value):
        return None
    token = str(value).strip().upper()
    if token.startswith("CHR"):
        token = token[3:]
    if token in CHR_ALIASES:
        return CHR_ALIASES[token]
    try:
        iv = int(token)
    except ValueError:
        return None
    if iv < 1:
        return None
    return iv


def preprocess_sumstats(df: pd.DataFrame, skip: float = 0.0) -> pd.DataFrame:
    """Convert CHR/POS/P and compute mlog10p."""
    out = df.copy()
    out["CHR"] = out["CHR"].map(_normalize_chr_token)
    out["POS"] = pd.to_numeric(out["POS"], errors="coerce")
    out["P"] = pd.to_numeric(out["P"], errors="coerce")
    out = out.dropna(subset=["CHR", "POS", "P"])
    out = out[(out["P"] > 0.0) & (out["P"] <= 1.0)]
    out["CHR"] = out["CHR"].astype(int)
    out["mlog10p"] = -np.log10(out["P"])
    if skip > 0:
        out = out[out["mlog10p"] >= float(skip)]
    return out.sort_values(["CHR", "POS"]).reset_index(drop=True)
