from pathlib import Path

from gwaspeek.interactive import (
    PAN_FRAC_COARSE,
    PAN_FRAC_FINE,
    ZOOM_IN_COARSE,
    ZOOM_IN_FINE,
    Viewport,
    _apply_key,
    _apply_mouse_wheel,
    _active_build,
    _parse_gtf_seqname_to_chrom,
    _parse_mouse_wheel,
    _plot_anchor_ratio,
    _protein_coding_track_in_view,
    _remap_viewport_to_build,
    parse_region,
    region_to_window,
)
from gwaspeek.io import load_sumstats
from gwaspeek.manhattan import genome_layout
from gwaspeek.plot_state import prepare_plot_dataset
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


def test_viewport_zoom_anchor_ratio_keeps_target_under_cursor() -> None:
    vp = Viewport(global_min=0.0, global_max=1_000_000.0, start=100_000.0, end=900_000.0, min_window_bp=1.0)
    anchor_ratio = 0.25
    anchor_before = vp.start + anchor_ratio * vp.width
    vp.zoom(2.0, anchor_ratio=anchor_ratio)
    anchor_after = vp.start + anchor_ratio * vp.width
    assert abs(anchor_after - anchor_before) < 1e-6


def test_apply_key_fine_pan_moves_less_than_coarse() -> None:
    _, _, x_min, x_max = _layout()
    vp_f = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max, min_window_bp=1.0)
    vp_f.zoom(10.0)
    vp_c = Viewport(
        global_min=x_min,
        global_max=x_max,
        start=vp_f.start,
        end=vp_f.end,
        min_window_bp=1.0,
    )
    w = vp_f.width
    start0 = vp_f.start
    _apply_key("d", vp_f, True, "37")
    _apply_key("D", vp_c, True, "37")
    assert vp_f.start == start0 + PAN_FRAC_FINE * w
    assert vp_c.start == start0 + PAN_FRAC_COARSE * w


def test_apply_key_fine_zoom_is_gentler_than_coarse() -> None:
    _, _, x_min, x_max = _layout()
    vp_f = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max, min_window_bp=1.0)
    vp_c = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max, min_window_bp=1.0)
    w0 = vp_f.width
    _apply_key("s", vp_f, True, "37")
    _apply_key("S", vp_c, True, "37")
    assert abs(vp_f.width - w0 / ZOOM_IN_FINE) < 1e-6
    assert abs(vp_c.width - w0 / ZOOM_IN_COARSE) < 1e-6


def test_apply_key_arrow_escape_does_not_pan_or_zoom() -> None:
    _, _, x_min, x_max = _layout()
    vp = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max, min_window_bp=1.0)
    vp.zoom(10.0)
    start, end, w = vp.start, vp.end, vp.width
    _apply_key("\x1b[C", vp, True, "37")
    assert vp.start == start and vp.end == end
    _apply_key("\x1b[A", vp, True, "37")
    assert vp.width == w


def test_apply_key_plus_minus_aliases_match_fine_zoom() -> None:
    _, _, x_min, x_max = _layout()
    zoom_in = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max, min_window_bp=1.0)
    zoom_out = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max, min_window_bp=1.0)
    w0 = zoom_in.width
    _apply_key("+", zoom_in, True, "37")
    _apply_key("-", zoom_out, True, "37")
    assert zoom_in.width < w0
    assert zoom_out.width == zoom_out.global_max - zoom_out.global_min


def test_parse_mouse_wheel_sgr_sequences() -> None:
    up = _parse_mouse_wheel("\x1b[<64;40;12M")
    down = _parse_mouse_wheel("\x1b[<65;41;13M")
    assert up is not None and up.direction == "up" and up.col == 40 and up.row == 12
    assert down is not None and down.direction == "down" and down.col == 41 and down.row == 13
    assert _parse_mouse_wheel("\x1b[<0;10;10M") is None


def test_plot_anchor_ratio_uses_plot_columns_only() -> None:
    ratio = _plot_anchor_ratio(col_1based=31, row_1based=10, frame_width=100, frame_height=28)
    assert ratio is not None
    assert abs(ratio - (23.0 / 92.0)) < 1e-6
    assert _plot_anchor_ratio(col_1based=3, row_1based=10, frame_width=100, frame_height=28) is None
    assert _plot_anchor_ratio(col_1based=31, row_1based=30, frame_width=100, frame_height=28) is None


def test_apply_mouse_wheel_zooms_at_cursor_anchor() -> None:
    vp = Viewport(global_min=0.0, global_max=1_000_000.0, start=100_000.0, end=900_000.0, min_window_bp=1.0)
    ratio = _plot_anchor_ratio(col_1based=31, row_1based=10, frame_width=100, frame_height=28)
    assert ratio is not None
    anchor_before = vp.start + ratio * vp.width
    applied = _apply_mouse_wheel("\x1b[<64;31;10M", vp, frame_width=100, frame_height=28)
    anchor_after = vp.start + ratio * vp.width
    assert applied is True
    assert vp.width < 800_000.0
    assert abs(anchor_after - anchor_before) < 1e-6


def test_apply_key_toggles_lead_annotation() -> None:
    _, _, x_min, x_max = _layout()
    vp = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max)
    keep, show_lead, track_mode = _apply_key("l", vp, True, "37")
    assert keep is True
    assert show_lead is False
    assert track_mode == "37"


def test_apply_key_non_human_ignores_gene_track_toggle() -> None:
    _, _, x_min, x_max = _layout()
    vp = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max)
    keep, show_lead, track_mode = _apply_key("t", vp, True, "off", non_human=True)
    assert keep is True
    assert track_mode == "off"


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


def test_active_build_uses_track_mode_or_base_build() -> None:
    assert _active_build("37", "38") == "37"
    assert _active_build("38", "37") == "38"
    assert _active_build("off", "37") == "37"
    assert _active_build("off", "38") == "38"


def test_remap_viewport_to_build_preserves_chr_positions() -> None:
    raw = load_sumstats(str(FIXTURE))
    clean = preprocess_sumstats(raw, skip=5.0)
    ds37 = prepare_plot_dataset(clean, build="37")
    ds38 = prepare_plot_dataset(clean, build="38")
    start37 = ds37.layout.offsets[1] + 1_000_000.0
    end37 = ds37.layout.offsets[2] + 2_000_000.0
    vp = Viewport(global_min=ds37.layout.x_min, global_max=ds37.layout.x_max, start=start37, end=end37, min_window_bp=1.0)
    _remap_viewport_to_build(vp, ds37, ds38)
    assert abs(vp.start - (ds38.layout.offsets[1] + 1_000_000.0)) < 1e-6
    assert abs(vp.end - (ds38.layout.offsets[2] + 2_000_000.0)) < 1e-6
    assert vp.global_min == ds38.layout.x_min
    assert vp.global_max == ds38.layout.x_max


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
