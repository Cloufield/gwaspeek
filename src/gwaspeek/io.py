from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd


def _normalize_sep(sep: str) -> str:
    # Shell usage sometimes passes escaped separators literally (e.g. "\\t").
    escape_map = {
        r"\t": "\t",
        r"\n": "\n",
        r"\r": "\r",
    }
    return escape_map.get(sep, sep)


@lru_cache(maxsize=1)
def _column_alias_table() -> Mapping[str, object]:
    path = resources.files("gwaspeek").joinpath("data/column_aliases_from_formatbook.json")
    return json.loads(path.read_text(encoding="utf-8"))


def _header_lower_index(headers: Sequence[str]) -> Dict[str, str]:
    """First occurrence wins for case-insensitive lookup."""
    out: Dict[str, str] = {}
    for h in headers:
        key = str(h).strip().lower()
        if key not in out:
            out[key] = str(h)
    return out


def _resolve_explicit(col: str, headers: Sequence[str], lower_index: Dict[str, str]) -> str:
    if col in headers:
        return col
    key = col.strip().lower()
    if key in lower_index:
        return lower_index[key]
    preview = ", ".join(repr(c) for c in list(headers)[:12])
    raise ValueError(f"Column {col!r} not found in input header (first columns: {preview}…)")


def _pick_alias(headers_set: set[str], lower_index: Dict[str, str], aliases: Sequence[str]) -> Optional[str]:
    """Prefer longer alias names first to avoid overly generic short tokens."""
    ranked = sorted(aliases, key=lambda a: (-len(str(a)), str(a).lower()))
    for alias in ranked:
        a = str(alias)
        if a in headers_set:
            return a
        key = a.strip().lower()
        if key in lower_index:
            return lower_index[key]
    return None


def detect_sumstat_columns(
    headers: Sequence[str],
    *,
    chrom_col: Optional[str] = None,
    pos_col: Optional[str] = None,
    p_col: Optional[str] = None,
    mlog10p_col: Optional[str] = None,
) -> Tuple[str, str, Optional[str], Optional[str]]:
    """
    Map file column names to roles CHR, POS, and either P or MLOG10P.

    Returns (chrom_src, pos_src, p_src, mlog10p_src) where exactly one of
    (p_src, mlog10p_src) is not None.
    """
    headers_list = [str(h) for h in headers]
    headers_set = set(headers_list)
    lower_index = _header_lower_index(headers_list)
    spec = _column_alias_table()
    raw_aliases = spec.get("aliases", {})
    aliases_chr = list(raw_aliases.get("CHR", ["CHR"]))
    aliases_pos = list(raw_aliases.get("POS", ["POS"]))
    aliases_p = list(raw_aliases.get("P", ["P"]))
    aliases_m = list(raw_aliases.get("MLOG10P", ["MLOG10P"]))

    if chrom_col is not None:
        chrom_src = _resolve_explicit(chrom_col, headers_list, lower_index)
    else:
        chrom_src = _pick_alias(headers_set, lower_index, aliases_chr)
        if chrom_src is None:
            raise ValueError("Could not auto-detect chromosome column (CHR). Set --chr.")

    if pos_col is not None:
        pos_src = _resolve_explicit(pos_col, headers_list, lower_index)
    else:
        pos_src = _pick_alias(headers_set, lower_index, aliases_pos)
        if pos_src is None:
            raise ValueError("Could not auto-detect position column (POS). Set --pos.")

    p_src: Optional[str] = None
    mlog_src: Optional[str] = None

    if p_col is not None:
        p_src = _resolve_explicit(p_col, headers_list, lower_index)
    if mlog10p_col is not None:
        mlog_src = _resolve_explicit(mlog10p_col, headers_list, lower_index)

    if p_src is None and mlog_src is None:
        p_src = _pick_alias(headers_set, lower_index, aliases_p)
        if p_src is None:
            mlog_src = _pick_alias(headers_set, lower_index, aliases_m)
        if p_src is None and mlog_src is None:
            raise ValueError(
                "Could not auto-detect P-value or -log10(P) column. "
                "Set --p and/or --mlog10p."
            )
    elif p_src is not None and mlog_src is not None:
        raise ValueError("Specify only one of --p or --mlog10p, not both.")

    return chrom_src, pos_src, p_src, mlog_src


def load_sumstats(
    path: str,
    sep: str = "\t",
    chrom_col: Optional[str] = None,
    pos_col: Optional[str] = None,
    p_col: Optional[str] = None,
    mlog10p_col: Optional[str] = None,
) -> pd.DataFrame:
    """Load GWAS summary statistics and normalize key column names to CHR, POS, and P or MLOG10P."""
    normalized_sep = _normalize_sep(sep)
    header_df = pd.read_csv(path, sep=normalized_sep, nrows=0, low_memory=False)
    headers = list(header_df.columns)

    chrom_src, pos_src, p_src, mlog_src = detect_sumstat_columns(
        headers,
        chrom_col=chrom_col,
        pos_col=pos_col,
        p_col=p_col,
        mlog10p_col=mlog10p_col,
    )

    usecols: List[str] = [chrom_src, pos_src]
    if p_src is not None:
        usecols.append(p_src)
    else:
        if mlog_src is None:
            raise ValueError("Internal error: MLOG10P column missing after detection.")
        usecols.append(mlog_src)

    required_set = set(usecols)
    df = pd.read_csv(
        path,
        sep=normalized_sep,
        low_memory=False,
        usecols=lambda c: c in required_set,
    )

    rename_map: Dict[str, str] = {
        chrom_src: "CHR",
        pos_src: "POS",
    }
    if p_src is not None:
        rename_map[p_src] = "P"
    else:
        rename_map[mlog_src] = "MLOG10P"

    return df.rename(columns=rename_map)
