from pathlib import Path

from gwaspeek.interactive import (
    Viewport,
    _apply_key,
    _parse_gtf_seqname_to_chrom,
    _protein_coding_track_in_view,
    parse_region,
    region_to_window,
)
from gwaspeek.io import load_sumstats
from gwaspeek.manhattan import genome_layout
from gwaspeek.preprocess import preprocess_sumstats


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"


def _layout():
    df = load_sumstats(str(FIXTURE))
    clean = preprocess_sumstats(df)
    _, offsets, chr_sizes, x_min, x_max = genome_layout(clean)
    return offsets, chr_sizes, x_min, x_max


def test_parse_region_chr_aliases() -> None:
    assert parse_region("chrX:10-20") == ("X", 10, 20)
    assert parse_region("M:15-30") == ("MT", 15, 30)


def test_region_to_window_maps_coordinates() -> None:
    offsets, chr_sizes, _, _ = _layout()
    start, end = region_to_window("1:100-200", offsets, chr_sizes)
    assert end > start
    assert start >= offsets[1]


def test_viewport_zoom_pan_and_reset() -> None:
    _, _, x_min, x_max = _layout()
    vp = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max, min_window_bp=1.0)
    original_width = vp.width
    vp.zoom(2.0)
    assert vp.width < original_width
    vp.pan(0.25 * vp.width)
    assert vp.end <= vp.global_max
    vp.reset()
    assert vp.start == x_min
    assert vp.end == x_max


def test_viewport_zoom_respects_default_10kb_floor() -> None:
    vp = Viewport(global_min=0.0, global_max=1_000_000.0, start=0.0, end=1_000_000.0)
    for _ in range(30):
        vp.zoom(2.0)
    assert vp.width >= 10_000.0


def test_apply_key_toggles_lead_annotation() -> None:
    _, _, x_min, x_max = _layout()
    vp = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max)
    keep, show_lead, track_mode = _apply_key("l", vp, True, "37")
    assert keep is True
    assert show_lead is False
    assert track_mode == "37"


def test_apply_key_cycles_gene_track_build_modes() -> None:
    _, _, x_min, x_max = _layout()
    vp = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max)
    keep, show_lead, track_mode = _apply_key("t", vp, True, "37")
    assert keep is True
    assert show_lead is True
    assert track_mode == "38"
    keep, show_lead, track_mode = _apply_key("t", vp, show_lead, track_mode)
    assert keep is True
    assert show_lead is True
    assert track_mode == "off"
    keep, show_lead, track_mode = _apply_key("t", vp, show_lead, track_mode)
    assert keep is True
    assert show_lead is True
    assert track_mode == "37"


def test_protein_coding_track_appears_for_small_single_chr_window() -> None:
    offsets, chr_sizes, _, _ = _layout()
    genes_by_chr = {1: [(100, 200, "GENE1"), (300000, 301000, "GENE2")]}
    start = offsets[1] + 50.0
    end = offsets[1] + 250.0
    track = _protein_coding_track_in_view(start, end, offsets, chr_sizes, genes_by_chr)
    assert len(track) == 1
    assert track[0][2] == "GENE1"


def test_protein_coding_track_when_viewport_crosses_chromosome_boundary() -> None:
    """Zoomed window can span cumulative chr boundary; still show genes on center chromosome."""
    offsets = {1: 0.0, 2: 1_000_000.0}
    chr_sizes = {1: 1_000_000, 2: 1_000_000}
    genes_by_chr = {
        1: [(999_400, 999_900, "TAIL1")],
        2: [(1000, 9000, "HEAD2")],
    }
    x_start = offsets[1] + 999_500.0
    x_end = offsets[2] + 5000.0
    track = _protein_coding_track_in_view(x_start, x_end, offsets, chr_sizes, genes_by_chr)
    # Midpoint lies on chr2; overlap with chr2 is [1e6, 1_005_000] -> bp 1..5000 includes HEAD2.
    names = {g[2] for g in track}
    assert "HEAD2" in names


def test_gtf_seqname_maps_nc_and_chr_prefixes() -> None:
    assert _parse_gtf_seqname_to_chrom("NC_000008.11") == 8
    assert _parse_gtf_seqname_to_chrom("chr8") == 8
    assert _parse_gtf_seqname_to_chrom("NC_012920.1") == 25
