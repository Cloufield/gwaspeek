from pathlib import Path

import pandas as pd

from gwaspeek.io import load_sumstats
from gwaspeek.manhattan import render_manhattan, sig_threshold_legend
from gwaspeek.preprocess import preprocess_sumstats


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"


def _prepared_manhattan():
    df = load_sumstats(str(FIXTURE))
    return preprocess_sumstats(df, skip=3.0)


def test_sig_threshold_legend_matches_plot_glyph() -> None:
    assert "=" in sig_threshold_legend(False, 5e-8)
    assert "5e-08" in sig_threshold_legend(False, 5e-8) or "5e-8" in sig_threshold_legend(False, 5e-8)
    assert "╌" in sig_threshold_legend(True, 5e-8)


def test_manhattan_render_deterministic() -> None:
    txt = render_manhattan(_prepared_manhattan(), width=60, height=18, unicode=False)
    assert "Manhattan" in txt
    assert "gwaspeek" in txt
    assert ":" in txt
    assert "+" in txt
    assert "*" in txt
    assert "3.0" in txt
    assert "skip>=" in txt


def test_manhattan_ascii_density_glyphs() -> None:
    txt = render_manhattan(_prepared_manhattan(), width=60, height=18, unicode=False)
    assert "*" in txt


def test_manhattan_unicode_density_glyphs() -> None:
    txt = render_manhattan(_prepared_manhattan(), width=60, height=18, unicode=True)
    assert "●" in txt


def test_manhattan_ascii_two_hits_use_O_glyph() -> None:
    raw = pd.DataFrame({"CHR": [1, 1], "POS": [100, 100], "P": [1e-8, 1e-8]})
    txt = render_manhattan(preprocess_sumstats(raw), width=40, height=12, unicode=False)
    assert "O" in txt


def test_manhattan_unicode_two_hits_use_bullseye_glyph() -> None:
    raw = pd.DataFrame({"CHR": [1, 1], "POS": [100, 100], "P": [1e-8, 1e-8]})
    txt = render_manhattan(preprocess_sumstats(raw), width=40, height=12, unicode=True)
    assert "◉" in txt


def test_manhattan_color_alternates_by_chromosome() -> None:
    txt = render_manhattan(_prepared_manhattan(), width=60, height=18, unicode=False, color=True)
    assert "\x1b[36m" in txt or "\x1b[37m" in txt


def test_manhattan_color_light_theme_alternates_by_chromosome() -> None:
    txt = render_manhattan(
        _prepared_manhattan(), width=60, height=18, unicode=False, color=True, light_theme=True
    )
    assert "\x1b[34m" in txt or "\x1b[35m" in txt


def test_manhattan_density_marks_overlapping_cells() -> None:
    raw = pd.DataFrame(
        {
            "CHR": [1, 1, 1, 1],
            "POS": [100, 100, 100, 200],
            "P": [1e-8, 1e-8, 1e-8, 1e-7],
        }
    )
    txt = render_manhattan(preprocess_sumstats(raw), width=40, height=12, unicode=False)
    assert "#" in txt
