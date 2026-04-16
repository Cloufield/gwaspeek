from pathlib import Path

from gwaspeek.io import load_sumstats
from gwaspeek.preprocess import preprocess_sumstats


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"


def test_preprocess_has_required_columns() -> None:
    df = load_sumstats(str(FIXTURE))
    out = preprocess_sumstats(df)
    assert {"CHR", "POS", "P", "mlog10p"}.issubset(out.columns)


def test_chr_normalization_ordering() -> None:
    df = load_sumstats(str(FIXTURE))
    out = preprocess_sumstats(df)
    assert out["CHR"].tolist() == sorted(out["CHR"].tolist())
    assert out["CHR"].max() == 25


def test_skip_filter() -> None:
    df = load_sumstats(str(FIXTURE))
    out = preprocess_sumstats(df, skip=3.0)
    assert (out["mlog10p"] >= 3.0).all()


def test_load_sumstats_accepts_escaped_tab_sep() -> None:
    df = load_sumstats(str(FIXTURE), sep=r"\t")
    assert {"CHR", "POS", "P"}.issubset(df.columns)
