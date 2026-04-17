from __future__ import annotations

import math
import os
import sys
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from gwaspeek.plot_state import (
    PlotDataset,
    build_genome_layout,
    prepare_plot_dataset,
    visible_mask,
)
from gwaspeek.terminal_canvas import CanvasStyle, TerminalCanvas
from gwaspeek.versioning import package_version


def _chr_token(chrom: int) -> str:
    return {23: "X", 24: "Y", 25: "MT"}.get(chrom, str(chrom))


def stdout_color_supported() -> bool:
    """True when stdout looks like a color-capable terminal (respects NO_COLOR)."""
    term = os.environ.get("TERM", "")
    return (
        sys.stdout.isatty()
        and term not in {"", "dumb"}
        and "NO_COLOR" not in os.environ
    )


def _chr_alternating_sgr(chrom: int, chr_order: Dict[int, int]) -> str:
    """ANSI SGR foreground for alternating colors by chromosome order (even/odd index)."""
    idx = chr_order.get(int(chrom))
    if idx is None:
        return "39"
    # Cyan (bluish) vs bright white; distinct for neighboring chromosomes.
    return "36" if (idx % 2) == 0 else "97"


def _ansi_paint_glyph(glyph: str, chrom: int, chr_order: Dict[int, int], color: bool) -> str:
    if not color:
        return glyph
    code = _chr_alternating_sgr(chrom, chr_order)
    return f"\x1b[{code}m{glyph}\x1b[0m"


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


def _format_pos_compact(pos: int) -> str:
    if pos >= 1_000_000:
        return f"{pos / 1e6:.1f}M"
    if pos >= 10_000:
        return f"{pos / 1e3:.0f}k"
    return str(pos)


def viewport_chr_label(
    x_lo: float,
    x_hi: float,
    offsets: Dict[int, float],
    chr_sizes: Dict[int, float],
) -> str:
    sorted_chroms = sorted(chr_sizes.keys())
    ca, pa = _cumulative_to_chr_pos(x_lo, sorted_chroms, offsets, chr_sizes)
    cb, pb = _cumulative_to_chr_pos(x_hi, sorted_chroms, offsets, chr_sizes)
    return f"{_chr_token(ca)}:{pa}-{_chr_token(cb)}:{pb}"


def _format_y_tick(v: float) -> str:
    return f"{float(v):.1f}"


def _nice_step(rough_step: float) -> float:
    if rough_step <= 0 or not math.isfinite(rough_step):
        return 0.1
    exp = math.floor(math.log10(rough_step))
    base = rough_step / (10**exp)
    if base <= 1.0:
        nice = 1.0
    elif base <= 2.0:
        nice = 2.0
    elif base <= 5.0:
        nice = 5.0
    else:
        nice = 10.0
    return nice * (10**exp)


def _dynamic_y_ticks(y_min: float, y_max: float, plot_height: int) -> List[float]:
    """Nice tick positions from y_min through a padded top; count scales with plot height."""
    if not math.isfinite(y_min) or not math.isfinite(y_max):
        return [float(y_min)] if math.isfinite(y_min) else [0.0]
    if y_max <= y_min:
        step = _nice_step(0.1)
        return [float(y_min), float(y_min) + step]
    span = y_max - y_min
    n_target = max(3, min(9, max(3, plot_height // 4)))
    step = _nice_step(span / max(2, n_target - 1))
    top = y_min + math.ceil((span * 1.02) / step) * step
    ticks: List[float] = []
    t = float(y_min)
    guard = 0
    while t <= top + 1e-9 and guard < 16:
        ticks.append(round(t, 8))
        t += step
        guard += 1
    return ticks if ticks else [float(y_min)]


def _chr_start_boundaries_in_view(
    x_min: float,
    x_max: float,
    sorted_chroms: List[int],
    offsets: Dict[int, float],
) -> List[Tuple[float, int]]:
    """Cumulative x where a new chromosome begins, strictly inside the visible x window."""
    out: List[Tuple[float, int]] = []
    for c in sorted_chroms[1:]:
        bx = float(offsets[c])
        if x_min < bx < x_max:
            out.append((bx, int(c)))
    return out


def _boundary_cx(
    bx: float,
    x_min: float,
    x_max: float,
    x_plot0: int,
    w: int,
) -> int:
    px = 0.0 if x_max <= x_min else (float(bx) - float(x_min)) / (float(x_max) - float(x_min))
    return x_plot0 + int(max(0.0, min(1.0, px)) * w)


def _draw_chr_transition_vlines(
    canvas: TerminalCanvas,
    boundaries: List[Tuple[float, int]],
    x_min: float,
    x_max: float,
    x_plot0: int,
    w: int,
    y_ranges: List[Tuple[int, int]],
    unicode: bool,
) -> None:
    v = "│" if unicode else "|"
    for bx, _ in boundaries:
        cx = _boundary_cx(bx, x_min, x_max, x_plot0, w)
        for y_lo, y_hi in y_ranges:
            y0 = max(0, y_lo)
            y1 = min(canvas.height - 1, y_hi)
            for y in range(y0, y1 + 1):
                canvas.set(cx, y, v)


def _mark_chr_boundaries_on_axis(
    canvas: TerminalCanvas,
    boundaries: List[Tuple[float, int]],
    x_min: float,
    x_max: float,
    x_plot0: int,
    w: int,
    axis_y: int,
    unicode: bool,
) -> None:
    """Mark chromosome boundaries on the horizontal x-axis (crosses the baseline)."""
    cross = "┼" if unicode else "+"
    for bx, _ in boundaries:
        cx = _boundary_cx(bx, x_min, x_max, x_plot0, w)
        if 0 <= axis_y < canvas.height:
            canvas.set(cx, axis_y, cross)


def _draw_dynamic_x_ticks(
    canvas: TerminalCanvas,
    x_min: float,
    x_max: float,
    global_x_min: float,
    global_x_max: float,
    sorted_chroms: List[int],
    offsets: Dict[int, float],
    chr_sizes: Dict[int, float],
    x_plot0: int,
    w: int,
    axis_y: int,
    label_y: int,
) -> None:
    full_span = global_x_max - global_x_min
    vis_span = x_max - x_min
    frac = (vis_span / full_span) if full_span > 0 else 1.0

    visible_chroms: List[int] = []
    for c in sorted_chroms:
        lo = offsets[c]
        hi = offsets[c] + chr_sizes[c]
        if hi >= x_min and lo <= x_max:
            visible_chroms.append(c)

    one_chrom = len(visible_chroms) == 1
    use_chr_centers = (not one_chrom) and frac >= 0.12

    ticks: List[Tuple[float, str]] = []
    if use_chr_centers:
        for c in visible_chroms:
            center = offsets[c] + (chr_sizes[c] / 2.0)
            if x_min <= center <= x_max:
                ticks.append((center, _chr_token(c)))
    if not ticks:
        n = max(2, min(7, max(1, w // 14)))
        for i in range(n):
            t = i / (n - 1) if n > 1 else 0.5
            x_c = x_min + t * (x_max - x_min)
            cc, pp = _cumulative_to_chr_pos(x_c, sorted_chroms, offsets, chr_sizes)
            ticks.append((x_c, f"{_chr_token(cc)}:{_format_pos_compact(pp)}"))

    used_cols: set[int] = set()
    for x_c, label in sorted(ticks, key=lambda z: z[0]):
        px = 0.0 if x_max <= x_min else (float(x_c) - x_min) / (x_max - x_min)
        tick_x = x_plot0 + int(max(0.0, min(1.0, px)) * w)
        if tick_x in used_cols:
            continue
        used_cols.add(tick_x)
        canvas.set(tick_x, axis_y, "+")
        start_x = max(x_plot0, tick_x - (len(label) // 2))
        for j, ch in enumerate(label):
            if start_x + j < canvas.width:
                canvas.set(start_x + j, label_y, ch)


def _draw_y_tick_labels(
    canvas: TerminalCanvas,
    ticks: List[float],
    y_min: float,
    y_scale: float,
    y0: int,
    h: int,
) -> None:
    """Draw -log10(P) tick labels in columns 0..5, right-aligned next to the y-axis."""
    denom = float(y_scale) - float(y_min)
    last_cy: int | None = None
    for yt in ticks:
        if denom <= 0:
            py_ratio = 0.0
        else:
            py_ratio = (float(yt) - float(y_min)) / denom
        cy = y0 - int(max(0.0, min(1.0, py_ratio)) * h)
        cy = max(1, min(y0, cy))
        if last_cy is not None and cy == last_cy:
            continue
        last_cy = cy
        label = _format_y_tick(yt).rjust(6)[:6]
        for i, ch in enumerate(label):
            canvas.set(i, cy, ch)


def _density_glyph(count: int, unicode: bool) -> str:
    if count <= 1:
        return "●" if unicode else "*"
    if count == 2:
        return "◉" if unicode else "O"
    return "█" if unicode else "#"


def density_legend(unicode: bool) -> str:
    if unicode:
        return "density 1x ●  2x ◉  3+x █"
    return "density 1x *  2x O  3+x #"


def _draw_gene_track(
    canvas: TerminalCanvas,
    genes: List[Tuple[float, float, str]],
    x_min: float,
    x_max: float,
    x_plot0: int,
    w: int,
    lane_rows: List[int],
    unicode: bool,
) -> None:
    if not genes or not lane_rows:
        return
    seg_ch = "━" if unicode else "-"
    left_cap = "├" if unicode else "|"
    right_cap = "┤" if unicode else "|"
    lane_pad = 1
    short_label_min_inner = 4

    # Greedy interval packing so overlapping genes use different lanes.
    lane_last_end: List[int] = [-10**9 for _ in lane_rows]
    label_end_by_row: Dict[int, int] = {}
    pixel_genes: List[Tuple[int, int, str]] = []
    for g_start, g_end, g_name in genes[:80]:
        if g_end < x_min or g_start > x_max:
            continue
        px0 = 0.0 if x_max <= x_min else (max(g_start, x_min) - x_min) / (x_max - x_min)
        px1 = 0.0 if x_max <= x_min else (min(g_end, x_max) - x_min) / (x_max - x_min)
        x0 = x_plot0 + int(max(0.0, min(1.0, px0)) * w)
        x1 = x_plot0 + int(max(0.0, min(1.0, px1)) * w)
        if x1 <= x0:
            x1 = min(canvas.width - 1, x0 + 1)
        pixel_genes.append((x0, x1, g_name))

    pixel_genes.sort(key=lambda g: (g[0], g[1]))
    for x0, x1, g_name in pixel_genes:
        lane_idx = 0
        while lane_idx < len(lane_rows) and x0 <= (lane_last_end[lane_idx] + lane_pad):
            lane_idx += 1
        if lane_idx >= len(lane_rows):
            continue
        lane_last_end[lane_idx] = x1
        y = lane_rows[lane_idx]
        canvas.set(x0, y, left_cap)
        canvas.set(x1, y, right_cap)
        for x in range(x0 + 1, x1):
            canvas.set(x, y, seg_ch)
        inner_start = x0 + 1
        inner_end = x1 - 1
        inner_width = inner_end - inner_start + 1
        if inner_width <= 0:
            continue
        if inner_width >= short_label_min_inner:
            label = g_name[: min(12, inner_width)]
            label_x = inner_start + max(0, (inner_width - len(label)) // 2)
            for i, ch in enumerate(label):
                x = label_x + i
                if inner_start <= x <= inner_end:
                    canvas.set(x, y, ch)
            continue

        # Short genes: place label on the row above/below to avoid tiny clipped inline names.
        label = g_name[:12]
        center = (x0 + x1) // 2
        label_x = max(x_plot0 + 1, center - (len(label) // 2))
        label_x = min(label_x, canvas.width - len(label) - 2)
        label_y = None
        for cand in (y - 1, y + 1):
            if cand < min(lane_rows) or cand > max(lane_rows):
                continue
            if label_x > label_end_by_row.get(cand, -10**9):
                label_y = cand
                break
        if label_y is None:
            continue
        for i, ch in enumerate(label):
            x = label_x + i
            if x_plot0 + 1 <= x <= canvas.width - 2:
                canvas.set(x, label_y, ch)
        label_end_by_row[label_y] = label_x + len(label) - 1


def _draw_gene_panel_frame(
    canvas: TerminalCanvas,
    x_plot0: int,
    x_right: int,
    panel_top: int,
    panel_bottom: int,
    unicode: bool,
) -> None:
    if panel_top >= panel_bottom:
        return
    h = "─" if unicode else "-"
    v = "│" if unicode else "|"
    tl = "┌" if unicode else "+"
    tr = "┐" if unicode else "+"
    bl = "└" if unicode else "+"
    br = "┘" if unicode else "+"
    for x in range(x_plot0, x_right + 1):
        canvas.set(x, panel_top, h)
        canvas.set(x, panel_bottom, h)
    for y in range(panel_top, panel_bottom + 1):
        canvas.set(x_plot0, y, v)
        canvas.set(x_right, y, v)
    canvas.set(x_plot0, panel_top, tl)
    canvas.set(x_right, panel_top, tr)
    canvas.set(x_plot0, panel_bottom, bl)
    canvas.set(x_right, panel_bottom, br)
    for i, ch in enumerate(" Genes "):
        pos = x_plot0 + 2 + i
        if pos < x_right:
            canvas.set(pos, panel_top, ch)


def _draw_lead_annotation(
    canvas: TerminalCanvas,
    lead_x: float,
    lead_y_ratio: float,
    lead_label: str,
    x_min: float,
    x_max: float,
    x_plot0: int,
    w: int,
    y0: int,
    h: int,
) -> None:
    px = 0.0 if x_max <= x_min else (lead_x - x_min) / (x_max - x_min)
    cx = x_plot0 + int(max(0.0, min(1.0, px)) * w)
    cy = y0 - int(max(0.0, min(1.0, lead_y_ratio)) * h)
    label = lead_label[:20]
    label_y = max(1, cy - 1)
    half = len(label) // 2
    start_x = min(max(x_plot0, cx - half), max(x_plot0, canvas.width - len(label) - 1))
    for i, ch in enumerate(label):
        if start_x + i < canvas.width:
            canvas.set(start_x + i, label_y, ch)


def _build_cumulative_x(df: pd.DataFrame) -> Tuple[pd.Series, Dict[int, float], Dict[int, float]]:
    layout = build_genome_layout(df)
    if len(df) == 0:
        xvals = pd.Series(dtype=float, index=df.index)
    else:
        xvals = pd.Series(
            df["POS"].to_numpy(dtype=float, copy=False) + df["CHR"].map(layout.offsets).to_numpy(dtype=float),
            index=df.index,
        )
    return xvals, layout.offsets, layout.chr_sizes


def genome_layout(df: pd.DataFrame) -> Tuple[pd.Series, Dict[int, float], Dict[int, float], float, float]:
    layout = build_genome_layout(df)
    if len(df) == 0:
        xvals = pd.Series(dtype=float, index=df.index)
    else:
        xvals = pd.Series(
            df["POS"].to_numpy(dtype=float, copy=False) + df["CHR"].map(layout.offsets).to_numpy(dtype=float),
            index=df.index,
        )
    return xvals, layout.offsets, layout.chr_sizes, layout.x_min, layout.x_max


def render_manhattan(
    df: pd.DataFrame,
    width: int = 100,
    height: int = 28,
    sig_level: float = 5e-8,
    ymax: float | None = None,
    x_start: float | None = None,
    x_end: float | None = None,
    title: str | None = None,
    unicode: bool = True,
    y_min: float = 5.0,
    lead_variant: Tuple[float, float, str] | None = None,
    gene_track: List[Tuple[float, float, str]] | None = None,
    force_gene_panel: bool = False,
    prepared: PlotDataset | None = None,
    visible_rows: np.ndarray | None = None,
    color: bool = False,
) -> str:
    data = prepared or prepare_plot_dataset(df)
    offsets = data.layout.offsets
    chr_sizes = data.layout.chr_sizes
    global_x_min = data.layout.x_min
    global_x_max = data.layout.x_max
    global_y_hi = float(np.nanmax(data.mlog10p)) if len(data.mlog10p) > 0 else float(y_min)
    sig_logp = float(-np.log10(sig_level))

    x_min = global_x_min if x_start is None else max(global_x_min, float(x_start))
    x_max = global_x_max if x_end is None else min(global_x_max, float(x_end))
    if x_max <= x_min:
        x_min = global_x_min
        x_max = global_x_max
    if title is None:
        title = f"Manhattan {viewport_chr_label(x_min, x_max, offsets, chr_sizes)}"
    mask = visible_rows if visible_rows is not None else visible_mask(data, x_min, x_max)
    if len(mask) != len(data.x):
        raise ValueError("Visible row mask does not match prepared dataset length.")

    if ymax is None:
        if np.any(mask):
            vis_hi = float(np.nanmax(data.mlog10p[mask]))
        else:
            vis_hi = global_y_hi
        y_floor = vis_hi
    else:
        y_floor = float(ymax)
    y_cap = max(y_floor, sig_logp, float(y_min) + 1e-9)

    canvas = TerminalCanvas(width=width, height=height, style=CanvasStyle(unicode=unicode))
    canvas.draw_axes()
    canvas.label_top_pair(title, f"gwaspeek v{package_version()}")

    x_plot0 = 7
    axis_y = canvas.height - 2
    label_y = canvas.height - 1
    has_gene_track = bool(force_gene_panel or gene_track)
    x_right = canvas.width - 1
    panel_top: int | None = None
    panel_bottom: int | None = None
    if has_gene_track and canvas.height >= 14:
        panel_bottom = axis_y - 1
        panel_top = panel_bottom - 3
    else:
        gene_track = []
        has_gene_track = False
    y0 = (panel_top - 2) if has_gene_track and panel_top is not None else (axis_y - 1)
    h = y0
    w = canvas.width - x_plot0 - 1
    y_ticks = _dynamic_y_ticks(float(y_min), y_cap, h)
    y_scale = max(y_ticks[-1], y_cap, sig_logp) if y_ticks else max(y_cap, sig_logp)
    denom = float(y_scale) - float(y_min)

    sorted_chroms = sorted(chr_sizes.keys())
    chr_order = {c: i for i, c in enumerate(sorted_chroms)}
    chr_boundaries = _chr_start_boundaries_in_view(x_min, x_max, sorted_chroms, offsets)

    cells: Dict[tuple[int, int], tuple[int, int, float]] = {}
    if np.any(mask):
        vis_x = data.x[mask]
        vis_y = data.mlog10p[mask]
        vis_chr = data.chrom[mask]
        vis_y_clipped = np.minimum(vis_y, y_scale)
        if x_max == x_min:
            px = np.zeros_like(vis_x)
        else:
            px = np.clip((vis_x - x_min) / (x_max - x_min), 0.0, 1.0)
        if denom <= 0:
            py = np.zeros_like(vis_y_clipped)
        else:
            py = np.clip((vis_y_clipped - float(y_min)) / denom, 0.0, 1.0)

        cx = x_plot0 + np.floor(px * w).astype(int)
        cy = y0 - np.floor(py * h).astype(int)
        flat = (cy * canvas.width) + cx
        order = np.lexsort((vis_y_clipped, flat))
        flat_sorted = flat[order]
        _, first_idx, counts = np.unique(flat_sorted, return_index=True, return_counts=True)
        best_idx = order[first_idx + counts - 1]
        for cell_flat, count, chrom in zip(flat[best_idx], counts, vis_chr[best_idx]):
            cell_y = int(cell_flat // canvas.width)
            cell_x = int(cell_flat % canvas.width)
            cells[(cell_x, cell_y)] = (int(count), int(chrom), 0.0)
    for (cx, cy), (count, chrom, _) in cells.items():
        glyph = _density_glyph(count, unicode)
        canvas.set(cx, cy, _ansi_paint_glyph(glyph, chrom, chr_order, color))

    if lead_variant is not None:
        lead_x, lead_y, lead_label = lead_variant
        if denom <= 0:
            lead_py = 0.0
        else:
            lead_py = max(0.0, min(1.0, (float(min(lead_y, y_scale)) - float(y_min)) / denom))
        _draw_lead_annotation(
            canvas=canvas,
            lead_x=float(lead_x),
            lead_y_ratio=lead_py,
            lead_label=lead_label,
            x_min=x_min,
            x_max=x_max,
            x_plot0=x_plot0,
            w=w,
            y0=y0,
            h=h,
        )

    # draw significance threshold line
    if denom <= 0:
        sig_ratio = 0.0
    else:
        sig_ratio = max(0.0, min(1.0, float((sig_logp - float(y_min)) / denom)))
    line_y = y0 - int(sig_ratio * h)
    for x in range(x_plot0, canvas.width - 1):
        canvas.set(x, line_y, "=" if not unicode else "╌")

    if has_gene_track and panel_top is not None and panel_bottom is not None:
        _draw_gene_panel_frame(canvas, x_plot0, x_right, panel_top, panel_bottom, unicode)
        lane_rows = list(range(panel_top + 1, panel_bottom))
        if gene_track:
            _draw_gene_track(
                canvas=canvas,
                genes=gene_track,
                x_min=x_min,
                x_max=x_max,
                x_plot0=x_plot0,
                w=w,
                lane_rows=lane_rows,
                unicode=unicode,
            )
        _draw_chr_transition_vlines(
            canvas,
            chr_boundaries,
            x_min,
            x_max,
            x_plot0,
            w,
            [(panel_top + 1, panel_bottom - 1)],
            unicode,
        )

    # draw chromosome ticks and labels at x-axis
    _draw_dynamic_x_ticks(
        canvas,
        x_min,
        x_max,
        global_x_min,
        global_x_max,
        sorted_chroms,
        offsets,
        chr_sizes,
        x_plot0,
        w,
        axis_y,
        label_y,
    )
    _mark_chr_boundaries_on_axis(
        canvas,
        chr_boundaries,
        x_min,
        x_max,
        x_plot0,
        w,
        axis_y,
        unicode,
    )

    _draw_y_tick_labels(canvas, y_ticks, float(y_min), y_scale, y0, h)
    return canvas.render()
