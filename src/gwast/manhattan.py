from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from gwast.terminal_canvas import CanvasStyle, TerminalCanvas


def _build_cumulative_x(df: pd.DataFrame) -> Tuple[pd.Series, Dict[int, float], Dict[int, float]]:
    chr_max = df.groupby("CHR")["POS"].max().sort_index()
    offset = 0.0
    offsets: Dict[int, float] = {}
    chr_sizes: Dict[int, float] = {}
    for chrom, max_pos in chr_max.items():
        ichrom = int(chrom)
        offsets[ichrom] = offset
        chr_sizes[ichrom] = float(max_pos)
        offset += float(max_pos)
    xvals = df.apply(lambda r: float(r["POS"]) + offsets[int(r["CHR"])], axis=1)
    return xvals, offsets, chr_sizes


def render_manhattan(
    df: pd.DataFrame,
    width: int = 100,
    height: int = 28,
    sig_level: float = 5e-8,
    ymax: float | None = None,
    unicode: bool = True,
) -> str:
    work = df.copy()
    work["x"], offsets, chr_sizes = _build_cumulative_x(work)
    x_min = float(work["x"].min())
    x_max = float(work["x"].max())
    y_vals = work["mlog10p"].to_numpy()
    y_cap = float(np.nanmax(y_vals) if ymax is None else ymax)
    y_cap = max(y_cap, -np.log10(sig_level))
    y_vals = np.clip(y_vals, 0, y_cap)

    points: List[Tuple[float, float]] = []
    for xv, yv in zip(work["x"], y_vals):
        px = 0.0 if x_max == x_min else (float(xv) - x_min) / (x_max - x_min)
        py = 0.0 if y_cap == 0 else float(yv) / y_cap
        points.append((px, py))

    canvas = TerminalCanvas(width=width, height=height, style=CanvasStyle(unicode=unicode))
    canvas.draw_axes()
    canvas.label_top("Manhattan plot")
    canvas.plot_points(points)

    # draw significance threshold line
    sig_y = max(0.0, min(1.0, float(-np.log10(sig_level) / y_cap)))
    y0 = canvas.height - 3
    h = y0
    line_y = y0 - int(sig_y * h)
    for x in range(7, canvas.width - 1):
        canvas.set(x, line_y, "=" if not unicode else "╌")

    # draw chromosome ticks and labels at x-axis
    axis_y = canvas.height - 2
    label_y = canvas.height - 1
    x0 = 7
    w = canvas.width - x0 - 1

    for chrom in sorted(chr_sizes):
        center = offsets[chrom] + (chr_sizes[chrom] / 2.0)
        px = 0.0 if x_max == x_min else (center - x_min) / (x_max - x_min)
        tick_x = x0 + int(max(0.0, min(1.0, px)) * w)
        canvas.set(tick_x, axis_y, "+")
        label = str(chrom if chrom <= 22 else {23: "X", 24: "Y", 25: "MT"}.get(chrom, chrom))
        start_x = max(x0, tick_x - (len(label) // 2))
        for i, ch in enumerate(label):
            if start_x + i < canvas.width:
                canvas.set(start_x + i, label_y, ch)
    return canvas.render()
