from __future__ import annotations

import gzip
import os
import re
import select
import shutil
import sys
import termios
import tty
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from gwaspeek.manhattan import density_legend, render_manhattan, sig_threshold_legend, viewport_chr_label
from gwaspeek.plot_state import PlotDataset, normalize_build, prepare_plot_dataset, visible_mask


REGION_RE = re.compile(r"^\s*([^:]+):(\d+)-(\d+)\s*$")
MOUSE_SGR_RE = re.compile(r"^\x1b\[<(\d+);(\d+);(\d+)([Mm])$")

# Pan/zoom: lowercase keys use fine steps; uppercase uses coarse (faster) steps.
PAN_FRAC_FINE = 0.10
PAN_FRAC_COARSE = 0.20
ZOOM_IN_FINE = 1.25
ZOOM_OUT_FINE = 1.0 / ZOOM_IN_FINE
ZOOM_IN_COARSE = 2.0
ZOOM_OUT_COARSE = 0.5


@dataclass
class Viewport:
    global_min: float
    global_max: float
    start: float
    end: float
    min_window_bp: float = 10_000.0

    @property
    def width(self) -> float:
        return self.end - self.start

    def _clamp(self) -> None:
        if self.start < self.global_min:
            shift = self.global_min - self.start
            self.start += shift
            self.end += shift
        if self.end > self.global_max:
            shift = self.end - self.global_max
            self.start -= shift
            self.end -= shift
        self.start = max(self.global_min, self.start)
        self.end = min(self.global_max, self.end)
        if self.end <= self.start:
            self.start = self.global_min
            self.end = self.global_max

    def zoom(self, factor: float, anchor_ratio: float = 0.5) -> None:
        if factor <= 0:
            raise ValueError("Zoom factor must be > 0")
        anchor_ratio = max(0.0, min(1.0, float(anchor_ratio)))
        anchor = self.start + anchor_ratio * self.width
        new_width = self.width / factor
        min_width = max(float(self.min_window_bp), 1.0)
        full_width = self.global_max - self.global_min
        new_width = min(max(new_width, min_width), full_width)
        self.start = anchor - anchor_ratio * new_width
        self.end = self.start + new_width
        self._clamp()

    def pan(self, delta: float) -> None:
        self.start += delta
        self.end += delta
        self._clamp()

    def set_window(self, start: float, end: float) -> None:
        if end <= start:
            raise ValueError("Region end must be greater than region start")
        self.start = start
        self.end = end
        self._clamp()

    def reset(self) -> None:
        self.start = self.global_min
        self.end = self.global_max


@dataclass(frozen=True)
class FrameSummary:
    region: str
    view_size: str
    n_vars: int
    n_chrs: int
    lead_label: str | None
    lead_gene: str | None
    gene_panel_active: bool


@dataclass(frozen=True)
class MouseWheelEvent:
    direction: str
    col: int
    row: int


@dataclass
class GeneTrackStore:
    gtf37_path: str
    gtf38_path: str
    _cache: Dict[str, Dict[int, List[Tuple[int, int, str]]]] = field(default_factory=dict)

    def path_for(self, mode: str) -> str:
        return self.gtf37_path if mode == "37" else self.gtf38_path

    def exists(self, mode: str) -> bool:
        return mode != "off" and os.path.exists(self.path_for(mode))

    def loaded(self, mode: str) -> bool:
        return mode in self._cache

    def get(self, mode: str) -> Dict[int, List[Tuple[int, int, str]]]:
        if mode == "off":
            return {}
        if mode not in self._cache:
            self._cache[mode] = _load_protein_coding_genes(self.path_for(mode))
        return self._cache[mode]


def parse_region(region: str) -> tuple[str, int, int]:
    match = REGION_RE.match(region)
    if not match:
        raise ValueError("Region must match chr:start-end (example: 3:100000-200000)")
    chrom_token = match.group(1).strip().upper()
    if chrom_token.startswith("CHR"):
        chrom_token = chrom_token[3:]
    start = int(match.group(2))
    end = int(match.group(3))
    if end <= start:
        raise ValueError("Region end must be greater than region start")
    if chrom_token == "M":
        chrom_token = "MT"
    return chrom_token, start, end


def region_to_window(region: str, offsets: dict[int, float], chr_sizes: dict[int, float]) -> tuple[float, float]:
    chrom_token, start, end = parse_region(region)
    alias = {"X": 23, "Y": 24, "MT": 25}
    try:
        chrom = alias.get(chrom_token, int(chrom_token))
    except ValueError as exc:
        raise ValueError(f"Invalid chromosome: {chrom_token}") from exc
    if chrom not in offsets or chrom not in chr_sizes:
        raise ValueError(f"Chromosome {chrom_token} not present in data")
    max_pos = int(chr_sizes[chrom])
    if start > max_pos:
        raise ValueError(f"Region start exceeds chromosome max position ({max_pos})")
    clamped_end = min(end, max_pos)
    if clamped_end <= start:
        raise ValueError("Region has no visible span after clamping to chromosome bounds")
    base = offsets[chrom]
    return base + float(start), base + float(clamped_end)


def _help_text() -> str:
    return (
        "Keys: q quit  h help  v vars  t 37|38|off  g jump  l lead  m theme  "
        "r reset  a/d pan (A/D fast)  w/s zoom (W/S fast)  wheel@cursor  +/-"
    )


def _chr_token(chrom: int) -> str:
    return {23: "X", 24: "Y", 25: "MT"}.get(chrom, str(chrom))


def _format_bp_span(span: float) -> str:
    bp = max(1, int(round(span)))
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.2f}Mb"
    if bp >= 1_000:
        return f"{bp / 1_000:.1f}kb"
    return f"{bp}bp"


def default_gtf_path() -> str:
    return str(resources.files("gwaspeek").joinpath("data/GRCh37_latest_genomic.gene_only.gtf.gz"))


def default_gtf38_path() -> str:
    return str(resources.files("gwaspeek").joinpath("data/GRCh38_latest_genomic.gene_only.gtf.gz"))


_NC_PREFIX_TO_CHROM: Dict[str, int] = {
    "NC_000001": 1,
    "NC_000002": 2,
    "NC_000003": 3,
    "NC_000004": 4,
    "NC_000005": 5,
    "NC_000006": 6,
    "NC_000007": 7,
    "NC_000008": 8,
    "NC_000009": 9,
    "NC_000010": 10,
    "NC_000011": 11,
    "NC_000012": 12,
    "NC_000013": 13,
    "NC_000014": 14,
    "NC_000015": 15,
    "NC_000016": 16,
    "NC_000017": 17,
    "NC_000018": 18,
    "NC_000019": 19,
    "NC_000020": 20,
    "NC_000021": 21,
    "NC_000022": 22,
    "NC_000023": 23,
    "NC_000024": 24,
    "NC_012920": 25,
}


def _parse_gtf_seqname_to_chrom(seqname: str) -> int | None:
    raw = seqname.strip()
    if not raw:
        return None
    token = raw.upper()
    if token.startswith("CHR"):
        token = token[3:]
    if token == "M":
        token = "MT"
    alias = {"X": 23, "Y": 24, "MT": 25}
    if token in alias:
        return alias[token]
    try:
        return int(token)
    except ValueError:
        pass
    if token.startswith("NC_") or token.startswith("NT_"):
        base = token.split(".", 1)[0]
        if base in _NC_PREFIX_TO_CHROM:
            return _NC_PREFIX_TO_CHROM[base]
    return None


def _cumulative_to_chr_pos(
    x_cum: float,
    sorted_chroms: List[int],
    offsets: Dict[int, float],
    chr_sizes: Dict[int, float],
) -> Tuple[int, int]:
    x_cum = float(x_cum)
    if not sorted_chroms:
        return 1, 1
    first = sorted_chroms[0]
    if x_cum <= offsets[first]:
        return first, 1
    last = sorted_chroms[-1]
    if x_cum >= offsets[last] + chr_sizes[last]:
        return last, max(1, int(chr_sizes[last]))
    for i, c in enumerate(sorted_chroms):
        lo = offsets[c]
        next_lo = offsets[sorted_chroms[i + 1]] if i + 1 < len(sorted_chroms) else offsets[c] + chr_sizes[c] + 1.0
        if lo <= x_cum < next_lo:
            pos = int(round(x_cum - lo))
            pos = max(1, min(int(chr_sizes[c]), pos))
            return c, pos
    return last, max(1, int(chr_sizes[last]))


def _parse_gtf_attrs(attr_field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for chunk in attr_field.split(";"):
        token = chunk.strip()
        if not token or " " not in token:
            continue
        key, raw = token.split(" ", 1)
        out[key] = raw.strip().strip('"')
    return out


@lru_cache(maxsize=2)
def _load_protein_coding_genes(gtf_path: str) -> Dict[int, List[Tuple[int, int, str]]]:
    genes: Dict[int, List[Tuple[int, int, str]]] = {}
    if not os.path.exists(gtf_path):
        return genes
    with gzip.open(gtf_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = _parse_gtf_seqname_to_chrom(parts[0])
            if chrom is None:
                continue
            attrs = _parse_gtf_attrs(parts[8])
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype")
            if gene_type != "protein_coding":
                continue
            start = int(parts[3])
            end = int(parts[4])
            name = attrs.get("gene_name") or attrs.get("gene_id") or "gene"
            genes.setdefault(chrom, []).append((start, end, name))
    for chrom in genes:
        genes[chrom].sort(key=lambda g: g[0])
    return genes


def _lead_variant_in_view(data: PlotDataset, mask: np.ndarray) -> Tuple[float, float, str] | None:
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return None
    lead_idx = int(idx[np.argmax(data.mlog10p[idx])])
    chrom = int(data.chrom[lead_idx])
    pos = int(data.pos[lead_idx])
    return float(data.x[lead_idx]), float(data.mlog10p[lead_idx]), f"lead {_chr_token(chrom)}:{pos}"


def _nearest_gene_label_for_x(
    x_cum: float,
    offsets: Dict[int, float],
    chr_sizes: Dict[int, float],
    genes_by_chr: Dict[int, List[Tuple[int, int, str]]],
) -> str | None:
    sorted_chroms = sorted(chr_sizes.keys())
    chrom, pos = _cumulative_to_chr_pos(x_cum, sorted_chroms, offsets, chr_sizes)
    genes = genes_by_chr.get(chrom, [])
    if not genes:
        return None
    best_name: str | None = None
    best_dist: int | None = None
    for g_start, g_end, g_name in genes:
        if g_start <= pos <= g_end:
            return g_name
        dist = g_start - pos if pos < g_start else pos - g_end
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_name = g_name
    return best_name


def _nearest_gene_for_chr_pos(
    chrom: int,
    pos: int,
    genes_by_chr: Dict[int, List[Tuple[int, int, str]]],
) -> str | None:
    genes = genes_by_chr.get(chrom, [])
    if not genes:
        return None
    best_name: str | None = None
    best_dist: int | None = None
    for g_start, g_end, g_name in genes:
        if g_start <= pos <= g_end:
            return g_name
        dist = g_start - pos if pos < g_start else pos - g_end
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_name = g_name
    return best_name


def _protein_coding_track_in_view(
    x_start: float,
    x_end: float,
    offsets: Dict[int, float],
    chr_sizes: Dict[int, float],
    genes_by_chr: Dict[int, List[Tuple[int, int, str]]],
) -> List[Tuple[float, float, str]]:
    sorted_chroms = sorted(chr_sizes.keys())
    c0, p0 = _cumulative_to_chr_pos(x_start, sorted_chroms, offsets, chr_sizes)
    c1, p1 = _cumulative_to_chr_pos(x_end, sorted_chroms, offsets, chr_sizes)

    def _genes_on_chrom(chrom: int, lo_bp: int, hi_bp: int) -> List[Tuple[float, float, str]]:
        span_bp = abs(hi_bp - lo_bp)
        if span_bp > 1_000_000:
            return []
        genes = genes_by_chr.get(chrom, [])
        if not genes:
            return []
        lo = min(lo_bp, hi_bp)
        hi = max(lo_bp, hi_bp)
        out: List[Tuple[float, float, str]] = []
        off = offsets[chrom]
        for g_start, g_end, g_name in genes:
            if g_end < lo or g_start > hi:
                continue
            out.append((off + float(g_start), off + float(g_end), g_name))
            if len(out) >= 25:
                break
        return out

    if c0 == c1:
        return _genes_on_chrom(c0, int(p0), int(p1))

    x_lo = float(min(x_start, x_end))
    x_hi = float(max(x_start, x_end))
    x_mid = 0.5 * (x_lo + x_hi)
    c_mid, _ = _cumulative_to_chr_pos(x_mid, sorted_chroms, offsets, chr_sizes)
    chr_lo = float(offsets[c_mid])
    chr_hi = chr_lo + float(chr_sizes[c_mid])
    vx0 = max(x_lo, chr_lo)
    vx1 = min(x_hi, chr_hi)
    if vx1 <= vx0:
        return []
    p_lo = int(round(vx0 - chr_lo))
    p_hi = int(round(vx1 - chr_lo))
    p_lo = max(1, p_lo)
    p_hi = min(int(chr_sizes[c_mid]), p_hi)
    return _genes_on_chrom(c_mid, p_lo, p_hi)


def _inspect_view(
    data: PlotDataset,
    viewport: Viewport,
    show_lead: bool,
    show_track: bool,
    genes_by_chr: Dict[int, List[Tuple[int, int, str]]],
    *,
    y_min: float,
    non_human: bool = False,
) -> tuple[np.ndarray, Tuple[float, float, str] | None, List[Tuple[float, float, str]], FrameSummary, str]:
    if non_human:
        show_track = False
    mask = visible_mask(data, viewport.start, viewport.end)
    idx = np.flatnonzero(mask)
    n_vars = int(idx.size)
    n_chrs = int(np.unique(data.chrom[idx]).size) if n_vars > 0 else 0
    lead_variant = _lead_variant_in_view(data, mask) if show_lead else None
    lead_gene = None
    if lead_variant is not None and genes_by_chr:
        lead_gene = _nearest_gene_label_for_x(lead_variant[0], data.layout.offsets, data.layout.chr_sizes, genes_by_chr)

    track_allowed = show_track and viewport.width <= 1_000_000.0
    gene_track = (
        _protein_coding_track_in_view(viewport.start, viewport.end, data.layout.offsets, data.layout.chr_sizes, genes_by_chr)
        if track_allowed and genes_by_chr
        else []
    )

    region = viewport_chr_label(viewport.start, viewport.end, data.layout.offsets, data.layout.chr_sizes)
    view_size = _format_bp_span(viewport.width)
    title = f"Manhattan {region} | view {view_size} | vars {n_vars} | skip>={float(y_min):.1f}"
    if n_chrs > 1:
        title = f"{title} | chrs {n_chrs}"
    if lead_variant is not None:
        if lead_gene:
            title = f"{title} [{lead_variant[2]} | gene {lead_gene}]"
        else:
            title = f"{title} [{lead_variant[2]}]"

    summary = FrameSummary(
        region=region,
        view_size=view_size,
        n_vars=n_vars,
        n_chrs=n_chrs,
        lead_label=lead_variant[2] if lead_variant is not None else None,
        lead_gene=lead_gene,
        gene_panel_active=track_allowed,
    )
    return mask, lead_variant, gene_track, summary, title


def _render_frame(
    data: PlotDataset,
    viewport: Viewport,
    show_lead: bool,
    show_track: bool,
    genes_by_chr: Dict[int, List[Tuple[int, int, str]]],
    width: int,
    height: int,
    sig_level: float,
    ymax: float | None,
    unicode: bool,
    y_min: float,
    color: bool = False,
    *,
    non_human: bool = False,
    light_theme: bool = False,
) -> tuple[str, FrameSummary]:
    mask, lead_variant, gene_track, summary, title = _inspect_view(
        data,
        viewport,
        show_lead,
        show_track,
        genes_by_chr,
        y_min=y_min,
        non_human=non_human,
    )
    frame = render_manhattan(
        data.df,
        width=width,
        height=height,
        sig_level=sig_level,
        ymax=ymax,
        x_start=viewport.start,
        x_end=viewport.end,
        title=title,
        unicode=unicode,
        y_min=y_min,
        lead_variant=lead_variant,
        gene_track=gene_track,
        force_gene_panel=summary.gene_panel_active,
        prepared=data,
        visible_rows=mask,
        color=color,
        light_theme=light_theme,
    )
    return frame, summary


def _clear_screen() -> None:
    sys.stdout.write("\x1b[H\x1b[2J")


def _enter_alt_screen() -> None:
    sys.stdout.write("\x1b[?1049h\x1b[?1000h\x1b[?1006h\x1b[?25l")


def _leave_alt_screen() -> None:
    sys.stdout.write("\x1b[?25h\x1b[?1006l\x1b[?1000l\x1b[?1049l")


def _parse_mouse_wheel(seq: str) -> MouseWheelEvent | None:
    match = MOUSE_SGR_RE.match(seq)
    if not match or match.group(4) != "M":
        return None
    code = int(match.group(1))
    if (code & 64) == 0:
        return None
    wheel_button = code & 0b11
    if wheel_button == 0:
        direction = "up"
    elif wheel_button == 1:
        direction = "down"
    else:
        return None
    return MouseWheelEvent(direction=direction, col=int(match.group(2)), row=int(match.group(3)))


def _plot_anchor_ratio(col_1based: int, row_1based: int, frame_width: int, frame_height: int) -> float | None:
    if row_1based < 1 or row_1based > frame_height:
        return None
    x_plot0 = 7
    w = frame_width - x_plot0 - 1
    if w <= 0:
        return None
    x = col_1based - 1
    if x < x_plot0 or x > x_plot0 + w:
        return None
    return (x - x_plot0) / w


def _apply_mouse_wheel(seq: str, viewport: Viewport, frame_width: int, frame_height: int) -> bool:
    event = _parse_mouse_wheel(seq)
    if event is None:
        return False
    anchor_ratio = _plot_anchor_ratio(event.col, event.row, frame_width, frame_height)
    if anchor_ratio is None:
        return False
    if event.direction == "up":
        viewport.zoom(ZOOM_IN_FINE, anchor_ratio=anchor_ratio)
    else:
        viewport.zoom(ZOOM_OUT_FINE, anchor_ratio=anchor_ratio)
    return True


def _read_key_nonblocking(timeout_sec: float = 0.1) -> str | None:
    ready, _, _ = select.select([sys.stdin], [], [], timeout_sec)
    if not ready:
        return None
    ch = sys.stdin.read(1)
    if ch != "\x1b":
        return ch
    seq = ch
    for _ in range(32):
        more_ready, _, _ = select.select([sys.stdin], [], [], 0.03)
        if not more_ready:
            break
        seq += sys.stdin.read(1)
        if seq.endswith(("A", "B", "C", "D")) or MOUSE_SGR_RE.match(seq):
            break
    return seq


def _next_track_mode(track_mode: str) -> str:
    if track_mode == "37":
        return "38"
    if track_mode == "38":
        return "off"
    return "37"


def _active_build(track_mode: str, base_build: str) -> str:
    return track_mode if track_mode in {"37", "38"} else base_build


def _remap_cumulative_position(x_cum: float, old_data: PlotDataset, new_data: PlotDataset) -> float:
    sorted_old = sorted(old_data.layout.chr_sizes.keys())
    chrom, pos = _cumulative_to_chr_pos(x_cum, sorted_old, old_data.layout.offsets, old_data.layout.chr_sizes)
    if chrom not in new_data.layout.offsets:
        return new_data.layout.x_min
    max_pos = int(new_data.layout.chr_sizes[chrom])
    return float(new_data.layout.offsets[chrom] + max(1, min(max_pos, pos)))


def _remap_viewport_to_build(viewport: Viewport, old_data: PlotDataset, new_data: PlotDataset) -> None:
    new_start = _remap_cumulative_position(viewport.start, old_data, new_data)
    new_end = _remap_cumulative_position(viewport.end, old_data, new_data)
    viewport.global_min = new_data.layout.x_min
    viewport.global_max = new_data.layout.x_max
    viewport.start = new_start
    viewport.end = new_end
    viewport._clamp()


def _apply_key(
    key: str,
    viewport: Viewport,
    show_lead: bool,
    track_mode: str,
    *,
    non_human: bool = False,
) -> Tuple[bool, bool, str]:
    if key in {"q", "Q"}:
        return False, show_lead, track_mode
    if key in {"r", "R"}:
        viewport.reset()
        return True, show_lead, track_mode
    if key in {"l", "L"}:
        return True, (not show_lead), track_mode
    if key in {"t", "T"}:
        if non_human:
            return True, show_lead, "off"
        return True, show_lead, _next_track_mode(track_mode)
    if key == "a":
        viewport.pan(-PAN_FRAC_FINE * viewport.width)
        return True, show_lead, track_mode
    if key == "A":
        viewport.pan(-PAN_FRAC_COARSE * viewport.width)
        return True, show_lead, track_mode
    if key == "d":
        viewport.pan(PAN_FRAC_FINE * viewport.width)
        return True, show_lead, track_mode
    if key == "D":
        viewport.pan(PAN_FRAC_COARSE * viewport.width)
        return True, show_lead, track_mode
    if key in {"w", "-", "_"}:
        viewport.zoom(ZOOM_OUT_FINE)
        return True, show_lead, track_mode
    if key == "W":
        viewport.zoom(ZOOM_OUT_COARSE)
        return True, show_lead, track_mode
    if key in {"s", "+", "="}:
        viewport.zoom(ZOOM_IN_FINE)
        return True, show_lead, track_mode
    if key == "S":
        viewport.zoom(ZOOM_IN_COARSE)
        return True, show_lead, track_mode
    return True, show_lead, track_mode


def _render_help_screen() -> str:
    lines = [
        "gwaspeek interactive help",
        "",
        "Navigation:",
        "  a/d          : pan (fine)",
        "  A/D          : pan (coarse)",
        "  w/s          : zoom out/in (fine)",
        "  W/S          : zoom out/in (coarse)",
        "  +/-          : zoom in/out (fine)",
        "  mouse wheel  : zoom at cursor position",
        "  g            : jump to region (chr:start-end)",
        "",
        "Toggles:",
        "  l : toggle lead variant annotation",
        "  t : cycle gene track mode (37 -> 38 -> off)",
        "  v : toggle variants-in-view list",
        "  h : toggle this help",
        "  m : toggle dark / light palette (plot + status colors)",
        "  r : reset to full genome view",
        "  q : quit",
        "",
        "Dense cells use heavier glyphs instead of dropping overlapping points.",
        "Press h (or v) again to return to plot.",
    ]
    return "\n".join(lines)


def _render_variants_view(
    data: PlotDataset,
    viewport: Viewport,
    mask: np.ndarray,
    summary: FrameSummary,
    genes_by_chr: Dict[int, List[Tuple[int, int, str]]],
    y_min: float,
    limit: int = 25,
) -> str:
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return f"Variants in view: none\nRegion: {summary.region}\nSkip threshold: -log10P >= {float(y_min):.2f}"
    order = idx[np.argsort(data.mlog10p[idx])[::-1][:limit]]
    lines = [
        f"Variants in view ({len(order)} shown, sorted by -log10P):",
        f"Region: {summary.region}",
        f"Skip threshold: -log10P >= {float(y_min):.2f}",
        "CHR\tPOS\tP\t-log10P\tnearest_gene",
    ]
    for i in order:
        chrom = int(data.chrom[i])
        pos = int(data.pos[i])
        nearest = _nearest_gene_for_chr_pos(chrom, pos, genes_by_chr) or "-"
        lines.append(f"{chrom}\t{pos}\t{float(data.p[i]):.3e}\t{float(data.mlog10p[i]):.3f}\t{nearest}")
    return "\n".join(lines)


def _genes_for_view(
    store: GeneTrackStore,
    track_mode: str,
    viewport: Viewport,
    view_mode: str,
    *,
    non_human: bool = False,
) -> Dict[int, List[Tuple[int, int, str]]]:
    if non_human or track_mode == "off":
        return {}
    if view_mode == "variants" or viewport.width <= 1_000_000.0:
        return store.get(track_mode)
    return {}


def _ansi_enabled(enabled: bool) -> bool:
    term = os.environ.get("TERM", "")
    return enabled and sys.stdout.isatty() and term not in {"", "dumb"} and "NO_COLOR" not in os.environ


def _ansi(text: str, code: str, enabled: bool) -> str:
    if not enabled:
        return text
    return f"\x1b[{code}m{text}\x1b[0m"


def _truncate_line(text: str, width: int) -> str:
    if width <= 0:
        return ""
    if len(text) <= width:
        return text
    if width <= 3:
        return text[:width]
    return f"{text[: width - 3].rstrip()}..."


def _format_track_state(
    track_mode: str,
    store: GeneTrackStore,
    summary: FrameSummary,
    *,
    non_human: bool = False,
) -> str:
    if non_human:
        return "track off (non-human)"
    if track_mode == "off":
        return "track off"
    if not store.exists(track_mode):
        return f"track {track_mode} missing"
    if summary.gene_panel_active:
        state = "loaded" if store.loaded(track_mode) else "loading"
        return f"track {track_mode} {state}"
    return f"track {track_mode} ready (zoom <=1Mb)"


def _status_line(
    summary: FrameSummary,
    build: str,
    track_mode: str,
    store: GeneTrackStore,
    notice: str | None,
    width: int,
    color: bool,
    *,
    light_theme: bool = False,
    non_human: bool = False,
) -> str:
    parts = [
        f"build {build}",
        f"region {summary.region}",
        f"view {summary.view_size}",
        f"vars {summary.n_vars}",
        f"chrs {summary.n_chrs}",
        _format_track_state(track_mode, store, summary, non_human=non_human),
    ]
    if summary.lead_label:
        lead = summary.lead_label
        if summary.lead_gene:
            lead = f"{lead} ({summary.lead_gene})"
        parts.append(lead)
    if notice:
        parts.append(notice)
    line = _truncate_line(" | ".join(parts), width)
    # Light palette: higher-contrast hues on pale backgrounds.
    if notice:
        sgr = "31" if light_theme else "33"
    else:
        sgr = "34" if light_theme else "36"
    return _ansi(line, sgr, color)


def _footer_help_line(
    width: int,
    unicode: bool,
    color: bool,
    sig_level: float,
    *,
    light_theme: bool = False,
) -> str:
    """Two-line footer; keys and density/sig legend share the same ANSI style."""
    line1 = _truncate_line(_help_text(), width)
    line2 = _truncate_line(
        f"{density_legend(unicode)}  |  {sig_threshold_legend(unicode, sig_level)}",
        width,
    )
    text = f"{line1}\n{line2}"
    # SGR 2 (faint) is often unsupported or invisible; use explicit greys (90 dark / 30 light).
    muted = "30" if light_theme else "90"
    return _ansi(text, muted, color)


def _prompt_for_region(
    fd: int,
    old_settings: list[object],
    viewport: Viewport,
    offsets: Dict[int, float],
    chr_sizes: Dict[int, float],
) -> str:
    sys.stdout.write("\x1b[?25h")
    sys.stdout.write("\x1b[H\x1b[2J")
    sys.stdout.write("Jump to region (chr:start-end, blank cancels): ")
    sys.stdout.flush()
    termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    try:
        region = sys.stdin.readline().strip()
    finally:
        tty.setcbreak(fd)
        sys.stdout.write("\x1b[?25l")
        sys.stdout.flush()
    if not region:
        return "Jump cancelled"
    try:
        start, end = region_to_window(region, offsets, chr_sizes)
        viewport.set_window(start, end)
    except ValueError as exc:
        return f"Jump failed: {exc}"
    return f"Jumped to {region}"


def run_interactive_manhattan(
    df: pd.DataFrame,
    width: int = 100,
    height: int = 28,
    sig_level: float = 5e-8,
    ymax: float | None = None,
    unicode: bool = True,
    y_min: float = 3.0,
    gtf_path: str | None = None,
    gtf38_path: str | None = None,
    build: str = "37",
    color: bool = True,
    non_human: bool = False,
) -> None:
    if not gtf_path:
        gtf_path = default_gtf_path()
    if not gtf38_path:
        gtf38_path = default_gtf38_path()

    base_build = normalize_build(build)
    datasets = {
        "37": prepare_plot_dataset(df, build="37", data_driven_lengths=non_human),
        "38": prepare_plot_dataset(df, build="38", data_driven_lengths=non_human),
    }
    show_lead = True
    track_mode = "off" if non_human else base_build
    view_mode = "plot"
    light_theme = False
    notice: str | None = None
    data = datasets[_active_build(track_mode, base_build)]
    offsets = data.layout.offsets
    chr_sizes = data.layout.chr_sizes
    viewport = Viewport(global_min=data.layout.x_min, global_max=data.layout.x_max, start=data.layout.x_min, end=data.layout.x_max)
    gene_store = GeneTrackStore(gtf37_path=gtf_path, gtf38_path=gtf38_path)

    if not sys.stdin.isatty():
        while True:
            active_build = _active_build(track_mode, base_build)
            data = datasets[active_build]
            offsets = data.layout.offsets
            chr_sizes = data.layout.chr_sizes
            genes_by_chr = _genes_for_view(
                gene_store, track_mode, viewport, view_mode, non_human=non_human
            )
            mask, _, _, summary, _ = _inspect_view(
                data,
                viewport,
                show_lead,
                track_mode != "off",
                genes_by_chr,
                y_min=y_min,
                non_human=non_human,
            )
            _clear_screen()
            if view_mode == "help":
                print(_render_help_screen())
            elif view_mode == "variants":
                print(_render_variants_view(data, viewport, mask, summary, genes_by_chr, y_min))
            else:
                frame, _ = _render_frame(
                    data,
                    viewport,
                    show_lead,
                    track_mode != "off",
                    genes_by_chr,
                    width,
                    height,
                    sig_level,
                    ymax,
                    unicode,
                    y_min,
                    color=False,
                    non_human=non_human,
                )
                print(frame)
            print(
                _status_line(
                    summary,
                    active_build,
                    track_mode,
                    gene_store,
                    notice,
                    width,
                    color=False,
                    light_theme=False,
                    non_human=non_human,
                )
            )
            print(
                _footer_help_line(width, unicode, color=False, sig_level=sig_level, light_theme=False)
            )
            key = sys.stdin.read(1)
            notice = None
            if not key:
                break
            if key in {"h", "H"}:
                view_mode = "plot" if view_mode == "help" else "help"
                continue
            if key in {"v", "V"}:
                view_mode = "plot" if view_mode == "variants" else "variants"
                continue
            old_build = active_build
            keep_going, show_lead, track_mode = _apply_key(
                key, viewport, show_lead, track_mode, non_human=non_human
            )
            new_build = _active_build(track_mode, base_build)
            if new_build != old_build:
                _remap_viewport_to_build(viewport, datasets[old_build], datasets[new_build])
                offsets = datasets[new_build].layout.offsets
                chr_sizes = datasets[new_build].layout.chr_sizes
            if not keep_going:
                break
        return

    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    color_enabled = _ansi_enabled(color)
    try:
        _enter_alt_screen()
        tty.setcbreak(fd)
        painted = False
        last_term: tuple[int, int] | None = None
        last_ui_snap: tuple[object, ...] | None = None
        while True:
            key = _read_key_nonblocking(timeout_sec=0.15)
            term_size = shutil.get_terminal_size(fallback=(width, height + 3))
            frame_width = max(20, term_size.columns)
            term_lines = term_size.lines
            # Leave room for status line plus two-line footer (no clipping / missing chrome).
            frame_height = max(8, term_lines - 3)
            term_sig = (frame_width, term_lines)
            resized = last_term is not None and term_sig != last_term
            last_term = term_sig

            # Apply input before paint so this frame matches the key (avoids one-key lag on m/h/v/...).
            active_build = _active_build(track_mode, base_build)
            if key is not None:
                offsets_pre = datasets[active_build].layout.offsets
                chr_sizes_pre = datasets[active_build].layout.chr_sizes
                if view_mode == "plot" and _apply_mouse_wheel(key, viewport, frame_width, frame_height):
                    pass
                elif key in {"h", "H"}:
                    view_mode = "plot" if view_mode == "help" else "help"
                elif key in {"v", "V"}:
                    view_mode = "plot" if view_mode == "variants" else "variants"
                elif key in {"g", "G"}:
                    view_mode = "plot"
                    notice = _prompt_for_region(fd, old_settings, viewport, offsets_pre, chr_sizes_pre)
                elif key in {"m", "M"}:
                    light_theme = not light_theme
                else:
                    old_build = active_build
                    keep_going, show_lead, track_mode = _apply_key(
                        key, viewport, show_lead, track_mode, non_human=non_human
                    )
                    new_build = _active_build(track_mode, base_build)
                    if new_build != old_build:
                        _remap_viewport_to_build(viewport, datasets[old_build], datasets[new_build])
                    if not keep_going:
                        break

            active_build = _active_build(track_mode, base_build)
            data = datasets[active_build]

            ui_snap = (
                view_mode,
                show_lead,
                track_mode,
                light_theme,
                viewport.start,
                viewport.end,
                viewport.width,
            )
            if painted and key is None and not resized and ui_snap == last_ui_snap:
                continue

            genes_by_chr = _genes_for_view(
                gene_store, track_mode, viewport, view_mode, non_human=non_human
            )
            mask, _, _, summary, _ = _inspect_view(
                data,
                viewport,
                show_lead,
                track_mode != "off",
                genes_by_chr,
                y_min=y_min,
                non_human=non_human,
            )

            _clear_screen()
            if view_mode == "help":
                body = _render_help_screen()
            elif view_mode == "variants":
                body = _render_variants_view(data, viewport, mask, summary, genes_by_chr, y_min)
            else:
                body, summary = _render_frame(
                    data,
                    viewport,
                    show_lead,
                    track_mode != "off",
                    genes_by_chr,
                    frame_width,
                    frame_height,
                    sig_level,
                    ymax,
                    unicode,
                    y_min,
                    color=color_enabled,
                    non_human=non_human,
                    light_theme=light_theme,
                )
            sys.stdout.write(body)
            if not body.endswith("\n"):
                sys.stdout.write("\n")
            sys.stdout.write(
                _status_line(
                    summary,
                    active_build,
                    track_mode,
                    gene_store,
                    notice,
                    frame_width,
                    color_enabled,
                    light_theme=light_theme,
                    non_human=non_human,
                )
            )
            sys.stdout.write("\n")
            sys.stdout.write(
                _footer_help_line(
                    frame_width,
                    unicode,
                    color_enabled,
                    sig_level,
                    light_theme=light_theme,
                )
            )
            sys.stdout.flush()

            painted = True
            notice = None
            last_ui_snap = ui_snap
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        _leave_alt_screen()
        sys.stdout.flush()
