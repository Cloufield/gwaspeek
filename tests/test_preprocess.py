from pathlib import Path

from gwast.io import load_sumstats
from gwast.preprocess import preprocess_sumstats


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
