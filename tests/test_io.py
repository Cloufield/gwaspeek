from pathlib import Path

import pandas as pd
import pytest

from gwaspeek.io import detect_sumstat_columns, load_sumstats
from gwaspeek.preprocess import preprocess_sumstats


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"
ALT = Path(__file__).parent / "fixtures" / "sumstats_formatbook_headers.tsv"
MLOG = Path(__file__).parent / "fixtures" / "sumstats_mlog10p.tsv"


def test_detect_columns_default_fixture_headers() -> None:
    h = pd.read_csv(FIXTURE, sep="\t", nrows=0, low_memory=False).columns.tolist()
    c, p, pv, ml = detect_sumstat_columns(h)
    assert c == "CHR" and p == "POS" and pv == "P" and ml is None


def test_detect_columns_formatbook_style_headers() -> None:
    h = pd.read_csv(ALT, sep="\t", nrows=0, low_memory=False).columns.tolist()
    c, p, pv, ml = detect_sumstat_columns(h)
    assert c == "#CHROM" and p == "GENPOS" and pv == "pval_nominal" and ml is None


def test_load_autodetect_alt_headers_roundtrip() -> None:
    df = load_sumstats(str(ALT))
    assert set(df.columns) == {"CHR", "POS", "P"}
    clean = preprocess_sumstats(df)
    assert len(clean) == 7


def test_load_mlog10p_only() -> None:
    df = load_sumstats(str(MLOG))
    assert set(df.columns) == {"CHR", "POS", "MLOG10P"}
    clean = preprocess_sumstats(df)
    assert "mlog10p" in clean.columns
    assert (clean["P"] > 0).all() and (clean["P"] <= 1).all()


def test_explicit_p_and_mlog10p_conflict() -> None:
    h = ["CHR", "POS", "P", "LOG10P"]
    with pytest.raises(ValueError, match="only one"):
        detect_sumstat_columns(h, p_col="P", mlog10p_col="LOG10P")
