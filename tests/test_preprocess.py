from pathlib import Path

import pandas as pd

from gwaspeek.io import load_sumstats
from gwaspeek.manhattan import genome_layout
from gwaspeek.plot_state import CANONICAL_CHR_LENGTHS, prepare_plot_dataset
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


def test_skip_filter_keeps_full_chromosome_layout_metadata() -> None:
    raw = pd.DataFrame(
        {
            "CHR": [1, 1, 2, 2],
            "POS": [100, 1_000, 100, 250_000_000],
            "P": [1e-8, 0.5, 1e-8, 0.5],
        }
    )
    clean_full = preprocess_sumstats(raw, skip=0.0)
    clean_hits = preprocess_sumstats(raw, skip=5.0)
    _, offsets_full, chr_sizes_full, _, _ = genome_layout(clean_full)
    _, offsets_hits, chr_sizes_hits, _, _ = genome_layout(clean_hits)
    assert offsets_hits == offsets_full
    assert chr_sizes_hits == chr_sizes_full


def test_prepare_plot_dataset_data_driven_uses_observed_max_not_canonical_x() -> None:
    raw = pd.DataFrame(
        {
            "CHR": [23, 23],
            "POS": [5_000_000, 8_000_000],
            "P": [1e-8, 1e-8],
        }
    )
    clean = preprocess_sumstats(raw, skip=5.0)
    ds = prepare_plot_dataset(clean, build="37", data_driven_lengths=True)
    assert ds.layout.chr_sizes[23] == 8_000_000.0
    assert ds.layout.chr_sizes[23] < CANONICAL_CHR_LENGTHS["37"][23]


def test_prepare_plot_dataset_uses_canonical_chr_lengths_for_build() -> None:
    raw = pd.DataFrame(
        {
            "CHR": [1, 2],
            "POS": [100, 100],
            "P": [1e-8, 1e-8],
        }
    )
    clean = preprocess_sumstats(raw, skip=5.0)
    ds37 = prepare_plot_dataset(clean, build="37")
    ds38 = prepare_plot_dataset(clean, build="38")
    assert ds37.layout.chr_sizes[1] == CANONICAL_CHR_LENGTHS["37"][1]
    assert ds37.layout.chr_sizes[2] == CANONICAL_CHR_LENGTHS["37"][2]
    assert ds38.layout.chr_sizes[1] == CANONICAL_CHR_LENGTHS["38"][1]
    assert ds38.layout.chr_sizes[2] == CANONICAL_CHR_LENGTHS["38"][2]
    assert ds37.layout.x_min == 0.0
    assert ds37.layout.x_max == CANONICAL_CHR_LENGTHS["37"][1] + CANONICAL_CHR_LENGTHS["37"][2]
