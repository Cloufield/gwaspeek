from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pandas as pd

from gwast.terminal_canvas import CanvasStyle, TerminalCanvas


def _qq_points(mlog10p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    obs = np.sort(mlog10p)[::-1]
    n = len(obs)
    exp = -np.log10(np.arange(1, n + 1) / (n + 1))
    return exp, obs


def render_qq(
    df: pd.DataFrame,
    width: int = 100,
    height: int = 28,
    ymax: float | None = None,
    unicode: bool = True,
) -> str:
    exp, obs = _qq_points(df["mlog10p"].to_numpy())
    top = float(max(np.nanmax(exp), np.nanmax(obs)))
    y_cap = top if ymax is None else float(ymax)
    y_cap = max(y_cap, 1.0)

    points: List[Tuple[float, float]] = []
    for xv, yv in zip(exp, obs):
        px = min(float(xv) / y_cap, 1.0)
        py = min(float(yv) / y_cap, 1.0)
        points.append((px, py))

    canvas = TerminalCanvas(width=width, height=height, style=CanvasStyle(unicode=unicode))
    canvas.draw_axes()
    canvas.label_top("QQ plot")
    canvas.plot_diag()
    canvas.plot_points(points)
    return canvas.render()
