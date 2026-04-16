from pathlib import Path

from gwaspeek.io import load_sumstats
from gwaspeek.manhattan import render_manhattan
from gwaspeek.preprocess import preprocess_sumstats


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"


def _prepared_manhattan():
    df = load_sumstats(str(FIXTURE))
    return preprocess_sumstats(df, skip=5.0)


def test_manhattan_render_deterministic() -> None:
    txt = render_manhattan(_prepared_manhattan(), width=60, height=18, unicode=False)
    assert "Manhattan" in txt
    assert "gwaspeek" in txt
    assert ":" in txt
    assert "+" in txt
    assert "*" in txt
    assert "o" in txt
    assert "5.0" in txt


def test_manhattan_ascii_uses_filled_and_hollow_dots() -> None:
    txt = render_manhattan(_prepared_manhattan(), width=60, height=18, unicode=False)
    assert "*" in txt
    assert "o" in txt


def test_manhattan_unicode_uses_filled_and_hollow_dots() -> None:
    txt = render_manhattan(_prepared_manhattan(), width=60, height=18, unicode=True)
    assert "●" in txt
    assert "○" in txt
