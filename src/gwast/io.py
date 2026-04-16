from __future__ import annotations

from typing import Dict

import pandas as pd


def load_sumstats(
    path: str,
    sep: str = "\t",
    chrom_col: str = "CHR",
    pos_col: str = "POS",
    p_col: str = "P",
) -> pd.DataFrame:
    """Load GWAS summary statistics and normalize key column names."""
    df = pd.read_csv(path, sep=sep, low_memory=False)
    required = [chrom_col, pos_col, p_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    rename_map: Dict[str, str] = {
        chrom_col: "CHR",
        pos_col: "POS",
        p_col: "P",
    }
    return df.rename(columns=rename_map)
