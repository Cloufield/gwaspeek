from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
import pandas as pd


DEFAULT_BUILD = "37"
CANONICAL_CHR_LENGTHS: dict[str, dict[int, float]] = {
    "37": {
        1: 249_250_621.0,
        2: 243_199_373.0,
        3: 198_022_430.0,
        4: 191_154_276.0,
        5: 180_915_260.0,
        6: 171_115_067.0,
        7: 159_138_663.0,
        8: 146_364_022.0,
        9: 141_213_431.0,
        10: 135_534_747.0,
        11: 135_006_516.0,
        12: 133_851_895.0,
        13: 115_169_878.0,
        14: 107_349_540.0,
        15: 102_531_392.0,
        16: 90_354_753.0,
        17: 81_195_210.0,
        18: 78_077_248.0,
        19: 59_128_983.0,
        20: 63_025_520.0,
        21: 48_129_895.0,
        22: 51_304_566.0,
        23: 155_270_560.0,
        24: 59_373_566.0,
        25: 16_569.0,
    },
    "38": {
        1: 248_956_422.0,
        2: 242_193_529.0,
        3: 198_295_559.0,
        4: 190_214_555.0,
        5: 181_538_259.0,
        6: 170_805_979.0,
        7: 159_345_973.0,
        8: 145_138_636.0,
        9: 138_394_717.0,
        10: 133_797_422.0,
        11: 135_086_622.0,
        12: 133_275_309.0,
        13: 114_364_328.0,
        14: 107_043_718.0,
        15: 101_991_189.0,
        16: 90_338_345.0,
        17: 83_257_441.0,
        18: 80_373_285.0,
        19: 58_617_616.0,
        20: 64_444_167.0,
        21: 46_709_983.0,
        22: 50_818_468.0,
        23: 156_040_895.0,
        24: 57_227_415.0,
        25: 16_569.0,
    },
}


@dataclass(frozen=True)
class GenomeLayout:
    build: str
    offsets: Dict[int, float]
    chr_sizes: Dict[int, float]
    x_min: float
    x_max: float


@dataclass(frozen=True)
class PlotDataset:
    df: pd.DataFrame
    layout: GenomeLayout
    x: np.ndarray
    chrom: np.ndarray
    pos: np.ndarray
    p: np.ndarray
    mlog10p: np.ndarray


def _coerce_chr_max_pos(raw_map: object) -> Dict[int, float] | None:
    if not isinstance(raw_map, dict):
        return None
    out: Dict[int, float] = {}
    for chrom, max_pos in raw_map.items():
        try:
            out[int(chrom)] = float(max_pos)
        except (TypeError, ValueError):
            continue
    return dict(sorted(out.items())) if out else None


def normalize_build(build: str | None) -> str:
    token = str(build or DEFAULT_BUILD).strip().upper()
    aliases = {
        "37": "37",
        "38": "38",
        "GRCH37": "37",
        "GRCH38": "38",
        "HG19": "37",
        "HG38": "38",
    }
    return aliases.get(token, DEFAULT_BUILD)


def build_genome_layout(
    df: pd.DataFrame,
    build: str | None = None,
    *,
    data_driven_lengths: bool = False,
) -> GenomeLayout:
    resolved_build = normalize_build(build or str(df.attrs.get("genome_build") or DEFAULT_BUILD))
    chr_max_pos = _coerce_chr_max_pos(df.attrs.get("chr_max_pos"))
    if chr_max_pos is None:
        chr_max = df.groupby("CHR")["POS"].max().sort_index()
        chr_max_pos = {int(chrom): float(max_pos) for chrom, max_pos in chr_max.items()}
    if not chr_max_pos:
        raise ValueError("No variants available to build genome layout.")

    canonical = CANONICAL_CHR_LENGTHS.get(resolved_build, {})
    offsets: Dict[int, float] = {}
    chr_sizes: Dict[int, float] = {}
    offset = 0.0
    for chrom, max_pos in chr_max_pos.items():
        offsets[chrom] = offset
        if data_driven_lengths:
            span = max(1.0, float(max_pos))
        else:
            span = float(canonical.get(chrom, max_pos))
        chr_sizes[chrom] = span
        offset += span

    sorted_chroms = sorted(chr_sizes)
    first_chrom = sorted_chroms[0]
    last_chrom = sorted_chroms[-1]
    x_min = float(offsets[first_chrom])
    x_max = float(offsets[last_chrom] + chr_sizes[last_chrom])

    if len(df) == 0:
        return GenomeLayout(build=resolved_build, offsets=offsets, chr_sizes=chr_sizes, x_min=x_min, x_max=x_max)

    return GenomeLayout(
        build=resolved_build,
        offsets=offsets,
        chr_sizes=chr_sizes,
        x_min=x_min,
        x_max=x_max,
    )


def prepare_plot_dataset(
    df: pd.DataFrame,
    layout: GenomeLayout | None = None,
    build: str | None = None,
    *,
    data_driven_lengths: bool = False,
) -> PlotDataset:
    resolved_layout = layout or build_genome_layout(
        df, build=build, data_driven_lengths=data_driven_lengths
    )
    if len(df) == 0:
        x = np.empty(0, dtype=float)
        chrom = np.empty(0, dtype=int)
        pos = np.empty(0, dtype=float)
        p = np.empty(0, dtype=float)
        mlog10p = np.empty(0, dtype=float)
    else:
        x = df["POS"].to_numpy(dtype=float, copy=False) + df["CHR"].map(resolved_layout.offsets).to_numpy(dtype=float)
        chrom = df["CHR"].to_numpy(dtype=int, copy=False)
        pos = df["POS"].to_numpy(dtype=float, copy=False)
        p = df["P"].to_numpy(dtype=float, copy=False)
        mlog10p = df["mlog10p"].to_numpy(dtype=float, copy=False)
    return PlotDataset(
        df=df,
        layout=resolved_layout,
        x=x,
        chrom=chrom,
        pos=pos,
        p=p,
        mlog10p=mlog10p,
    )


def visible_mask(dataset: PlotDataset, x_start: float, x_end: float) -> np.ndarray:
    if len(dataset.x) == 0:
        return np.zeros(0, dtype=bool)
    return (dataset.x >= float(x_start)) & (dataset.x <= float(x_end))
