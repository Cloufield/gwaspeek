from __future__ import annotations

import gzip
import os
import re
import select
import shutil
import sys
import termios
import tty
from importlib import resources
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Tuple

import pandas as pd

from gwaspeek.manhattan import genome_layout, render_manhattan, viewport_chr_label


REGION_RE = re.compile(r"^\s*([^:]+):(\d+)-(\d+)\s*$")


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

    def zoom(self, factor: float) -> None:
        if factor <= 0:
            raise ValueError("Zoom factor must be > 0")
        center = (self.start + self.end) / 2.0
        new_width = self.width / factor
        min_width = max(float(self.min_window_bp), 1.0)
        full_width = self.global_max - self.global_min
        new_width = min(max(new_width, min_width), full_width)
        self.start = center - (new_width / 2.0)
        self.end = center + (new_width / 2.0)
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
    return "Keys: A/D pan  W/S zoom  l lead  t 37|38|off  h help  v vars  r reset  q quit"


def _status_line(track_mode: str, gtf37_path: str, gtf38_path: str) -> str:
    gtf37_state = "ok" if os.path.exists(gtf37_path) else "missing"
    gtf38_state = "ok" if os.path.exists(gtf38_path) else "missing"
    return f"Gene track: {track_mode} | GTF37: {gtf37_state} | GTF38: {gtf38_state}"


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
    # Primary nuclear chromosomes (GRCh37/GRCh38 RefSeq naming; version suffix stripped).
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
    "NC_000023": 23,  # X
    "NC_000024": 24,  # Y
    # Mitochondrion (rCRS; common in NCBI RefSeq human assemblies)
    "NC_012920": 25,
}


def _parse_gtf_seqname_to_chrom(seqname: str) -> int | None:
    """Map GTF seq_id to internal numeric chromosome key used by gwaspeek."""
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


def _lead_variant_in_view(
    df: pd.DataFrame,
    x_start: float,
    x_end: float,
    offsets: Dict[int, float],
) -> Tuple[float, float, str] | None:
    xvals = df["POS"].astype(float) + df["CHR"].map(offsets).astype(float)
    vis = df[(xvals >= x_start) & (xvals <= x_end)]
    if len(vis) == 0:
        return None
    idx = vis["mlog10p"].idxmax()
    row = df.loc[idx]
    lead_x = float(row["POS"]) + float(offsets[int(row["CHR"])])
    lead_y = float(row["mlog10p"])
    lead_label = f"lead {_chr_token(int(row['CHR']))}:{int(row['POS'])}"
    return lead_x, lead_y, lead_label


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

    # Viewport crosses a chromosome boundary in cumulative space (common when zoomed
    # near a telomere). Single-chr logic returns nothing while the UI still reserves
    # the gene panel — use the chromosome at the viewport center and intersect [x_start, x_end]
    # with that chromosome's cumulative span.
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


def _render_frame(
    df: pd.DataFrame,
    viewport: Viewport,
    offsets: Dict[int, float],
    chr_sizes: Dict[int, float],
    show_lead: bool,
    show_track: bool,
    genes_by_chr: Dict[int, List[Tuple[int, int, str]]],
    width: int,
    height: int,
    sig_level: float,
    ymax: float | None,
    unicode: bool,
    y_min: float,
) -> str:
    xvals = df["POS"].astype(float) + df["CHR"].map(offsets).astype(float)
    vis = df[(xvals >= viewport.start) & (xvals <= viewport.end)]
    n_vars = len(vis)
    n_chrs = int(vis["CHR"].nunique()) if n_vars > 0 else 0

    lead_variant = _lead_variant_in_view(df, viewport.start, viewport.end, offsets) if show_lead else None
    view_size = _format_bp_span(viewport.width)
    title = (
        f"Manhattan {viewport_chr_label(viewport.start, viewport.end, offsets, chr_sizes)}"
        f" | view {view_size} | vars {n_vars}"
    )
    if n_chrs > 1:
        title = f"{title} | chrs {n_chrs}"
    if show_lead and lead_variant is not None:
        nearest_gene = _nearest_gene_label_for_x(lead_variant[0], offsets, chr_sizes, genes_by_chr)
        if nearest_gene:
            title = f"{title} [{lead_variant[2]} | gene {nearest_gene}]"
        else:
            title = f"{title} [{lead_variant[2]}]"
    track_allowed = show_track and viewport.width <= 1_000_000.0
    gene_track = (
        _protein_coding_track_in_view(viewport.start, viewport.end, offsets, chr_sizes, genes_by_chr)
        if track_allowed
        else []
    )
    return render_manhattan(
        df,
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
        force_gene_panel=track_allowed,
    )


def _clear_screen() -> None:
    print("\x1bc", end="")


def _read_key_nonblocking(timeout_sec: float = 0.1) -> str | None:
    ready, _, _ = select.select([sys.stdin], [], [], timeout_sec)
    if not ready:
        return None
    ch = sys.stdin.read(1)
    if ch != "\x1b":
        return ch
    # Escape-prefixed sequences (e.g. unused keys); allow a short grace period
    # so a single quick key press is not dropped when bytes arrive staggered.
    seq = ch
    for _ in range(5):
        more_ready, _, _ = select.select([sys.stdin], [], [], 0.03)
        if not more_ready:
            break
        seq += sys.stdin.read(1)
        if seq.endswith(("A", "B", "C", "D")):
            break
    return seq


def _next_track_mode(track_mode: str) -> str:
    if track_mode == "37":
        return "38"
    if track_mode == "38":
        return "off"
    return "37"


def _apply_key(key: str, viewport: Viewport, show_lead: bool, track_mode: str) -> Tuple[bool, bool, str]:
    if key in {"q", "Q"}:
        return False, show_lead, track_mode
    if key in {"r", "R"}:
        viewport.reset()
        return True, show_lead, track_mode
    if key in {"l", "L"}:
        return True, (not show_lead), track_mode
    if key in {"t", "T"}:
        return True, show_lead, _next_track_mode(track_mode)
    if key in {"a", "A"}:
        viewport.pan(-0.20 * viewport.width)
        return True, show_lead, track_mode
    if key in {"d", "D"}:
        viewport.pan(0.20 * viewport.width)
        return True, show_lead, track_mode
    if key in {"w", "W"}:
        viewport.zoom(0.5)
        return True, show_lead, track_mode
    if key in {"s", "S"}:
        viewport.zoom(2.0)
        return True, show_lead, track_mode
    return True, show_lead, track_mode


def _render_help_screen() -> str:
    lines = [
        "gwaspeek interactive help",
        "",
        "Navigation:",
        "  A/D : pan",
        "  W/S : zoom out/in",
        "",
        "Toggles:",
        "  l : toggle lead variant annotation",
        "  t : cycle gene track mode (37 -> 38 -> off)",
        "  v : toggle variants-in-view list",
        "  h : toggle this help",
        "  r : reset to full genome view",
        "  q : quit",
        "",
        "Press h (or v) again to return to plot.",
    ]
    return "\n".join(lines)


def _render_variants_view(
    df: pd.DataFrame,
    viewport: Viewport,
    offsets: Dict[int, float],
    chr_sizes: Dict[int, float],
    genes_by_chr: Dict[int, List[Tuple[int, int, str]]],
    y_min: float,
    limit: int = 25,
) -> str:
    xvals = df["POS"].astype(float) + df["CHR"].map(offsets).astype(float)
    vis = df[(xvals >= viewport.start) & (xvals <= viewport.end)].copy()
    region = viewport_chr_label(viewport.start, viewport.end, offsets, chr_sizes)
    if len(vis) == 0:
        return f"Variants in view: none\nRegion: {region}\nSkip threshold: -log10P >= {float(y_min):.2f}"
    vis = vis.sort_values("mlog10p", ascending=False).head(limit)
    lines = [
        f"Variants in view ({len(vis)} shown, sorted by -log10P):",
        f"Region: {region}",
        f"Skip threshold: -log10P >= {float(y_min):.2f}",
        "CHR\tPOS\tP\t-log10P\tnearest_gene",
    ]
    for row in vis.itertuples(index=False):
        chrom = int(getattr(row, "CHR"))
        pos = int(getattr(row, "POS"))
        nearest = _nearest_gene_for_chr_pos(chrom, pos, genes_by_chr) or "-"
        p_val = float(getattr(row, "P"))
        mlog10p = float(getattr(row, "mlog10p"))
        lines.append(f"{chrom}\t{pos}\t{p_val:.3e}\t{mlog10p:.3f}\t{nearest}")
    return "\n".join(lines)


def run_interactive_manhattan(
    df: pd.DataFrame,
    width: int = 100,
    height: int = 28,
    sig_level: float = 5e-8,
    ymax: float | None = None,
    unicode: bool = True,
    y_min: float = 5.0,
    gtf_path: str | None = None,
    gtf38_path: str | None = None,
) -> None:
    if not gtf_path:
        gtf_path = default_gtf_path()
    if not gtf38_path:
        gtf38_path = default_gtf38_path()
    _, offsets, chr_sizes, x_min, x_max = genome_layout(df)
    genes_by_chr_37 = _load_protein_coding_genes(gtf_path)
    genes_by_chr_38 = _load_protein_coding_genes(gtf38_path)
    show_lead = True
    track_mode = "37"
    view_mode = "plot"
    viewport = Viewport(global_min=x_min, global_max=x_max, start=x_min, end=x_max)
    if not sys.stdin.isatty():
        # Non-TTY fallback (e.g. tests): consume one-char commands from stdin.
        while True:
            genes_by_chr = genes_by_chr_37 if track_mode == "37" else genes_by_chr_38 if track_mode == "38" else {}
            _clear_screen()
            if view_mode == "help":
                print(_render_help_screen())
            elif view_mode == "variants":
                print(_render_variants_view(df, viewport, offsets, chr_sizes, genes_by_chr, y_min))
            else:
                print(
                    _render_frame(
                        df,
                        viewport,
                        offsets,
                        chr_sizes,
                        show_lead,
                        track_mode != "off",
                        genes_by_chr,
                        width,
                        height,
                        sig_level,
                        ymax,
                        unicode,
                        y_min,
                    )
                )
            print()
            print(_help_text())
            print(_status_line(track_mode, gtf_path, gtf38_path))
            key = sys.stdin.read(1)
            if not key:
                break
            if key in {"h", "H"}:
                view_mode = "plot" if view_mode == "help" else "help"
                continue
            if key in {"v", "V"}:
                view_mode = "plot" if view_mode == "variants" else "variants"
                continue
            keep_going, show_lead, track_mode = _apply_key(key, viewport, show_lead, track_mode)
            if not keep_going:
                break
        return

    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setcbreak(fd)
        while True:
            genes_by_chr = genes_by_chr_37 if track_mode == "37" else genes_by_chr_38 if track_mode == "38" else {}
            term_size = shutil.get_terminal_size(fallback=(width, height + 4))
            frame_width = max(20, term_size.columns)
            frame_height = max(8, term_size.lines - 3)
            _clear_screen()
            if view_mode == "help":
                print(_render_help_screen())
            elif view_mode == "variants":
                print(_render_variants_view(df, viewport, offsets, chr_sizes, genes_by_chr, y_min))
            else:
                print(
                    _render_frame(
                        df,
                        viewport,
                        offsets,
                        chr_sizes,
                        show_lead,
                        track_mode != "off",
                        genes_by_chr,
                        frame_width,
                        frame_height,
                        sig_level,
                        ymax,
                        unicode,
                        y_min,
                    )
                )
            print(_help_text())
            print(_status_line(track_mode, gtf_path, gtf38_path))
            key = _read_key_nonblocking(timeout_sec=0.15)
            if key is None:
                continue
            if key in {"h", "H"}:
                view_mode = "plot" if view_mode == "help" else "help"
                continue
            if key in {"v", "V"}:
                view_mode = "plot" if view_mode == "variants" else "variants"
                continue
            keep_going, show_lead, track_mode = _apply_key(key, viewport, show_lead, track_mode)
            if not keep_going:
                break
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
