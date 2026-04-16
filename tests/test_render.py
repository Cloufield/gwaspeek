from pathlib import Path

from gwast.io import load_sumstats
from gwast.manhattan import render_manhattan
from gwast.preprocess import preprocess_sumstats
from gwast.qq import render_qq


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"


def _prepared():
    df = load_sumstats(str(FIXTURE))
    return preprocess_sumstats(df)


def test_manhattan_render_deterministic() -> None:
    txt = render_manhattan(_prepared(), width=60, height=18, unicode=False)
    assert "Manhattan plot" in txt
    assert "+" in txt
    assert "o" in txt


def test_qq_render_deterministic() -> None:
    txt = render_qq(_prepared(), width=60, height=18, unicode=False)
    assert "QQ plot" in txt
    assert "/" in txt
    assert "o" in txt
