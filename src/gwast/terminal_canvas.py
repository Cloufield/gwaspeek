from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Tuple


@dataclass
class CanvasStyle:
    unicode: bool = True

    @property
    def point(self) -> str:
        return "●" if self.unicode else "o"

    @property
    def axis_h(self) -> str:
        return "─" if self.unicode else "-"

    @property
    def axis_v(self) -> str:
        return "│" if self.unicode else "|"

    @property
    def diag(self) -> str:
        return "╱" if self.unicode else "/"


class TerminalCanvas:
    def __init__(self, width: int, height: int, style: CanvasStyle):
        if width < 20 or height < 8:
            raise ValueError("Canvas too small, use width>=20 and height>=8")
        self.width = width
        self.height = height
        self.style = style
        self.grid: List[List[str]] = [[" " for _ in range(width)] for _ in range(height)]

    def set(self, x: int, y: int, ch: str) -> None:
        if 0 <= x < self.width and 0 <= y < self.height:
            self.grid[y][x] = ch

    def draw_axes(self) -> None:
        x0 = 6
        y0 = self.height - 2
        for x in range(x0, self.width):
            self.set(x, y0, self.style.axis_h)
        for y in range(0, y0 + 1):
            self.set(x0, y, self.style.axis_v)
        self.set(x0, y0, "+")

    def plot_points(self, points: Iterable[Tuple[float, float]]) -> None:
        x0 = 7
        y0 = self.height - 3
        w = self.width - x0 - 1
        h = y0
        occupied = set()
        for px, py in points:
            cx = x0 + int(max(0.0, min(1.0, px)) * w)
            cy = y0 - int(max(0.0, min(1.0, py)) * h)
            if (cx, cy) in occupied:
                continue
            occupied.add((cx, cy))
            self.set(cx, cy, self.style.point)

    def plot_diag(self) -> None:
        x0 = 7
        y0 = self.height - 3
        w = self.width - x0 - 1
        h = y0
        steps = min(w, h)
        for i in range(steps + 1):
            x = x0 + i
            y = y0 - i
            self.set(x, y, self.style.diag)

    def label_top(self, text: str) -> None:
        label = text[: self.width]
        for i, ch in enumerate(label):
            self.set(i, 0, ch)

    def render(self) -> str:
        return "\n".join("".join(row).rstrip() for row in self.grid)
