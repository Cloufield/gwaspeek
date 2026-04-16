from __future__ import annotations

import argparse
import sys
from importlib.metadata import PackageNotFoundError, version

from gwaspeek.interactive import default_gtf_path, run_interactive_manhattan
from gwaspeek.io import load_sumstats
from gwaspeek.manhattan import render_manhattan
from gwaspeek.preprocess import preprocess_sumstats


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gwaspeek",
        description="gwaspeek: interactive terminal viewer for GWAS summary statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  gwaspeek tests/fixtures/sumstats_small.tsv\n"
            "  gwaspeek -i tests/fixtures/sumstats_small.tsv\n"
            "  gwaspeek -s tests/fixtures/sumstats_small.tsv\n"
        ),
    )
    parser.add_argument(
        "sumstats",
        nargs="?",
        metavar="FILE",
        help="Summary statistics (interactive Manhattan; same as -i FILE; default mode)",
    )
    parser.add_argument(
        "-s",
        "--static",
        dest="static_file",
        metavar="FILE",
        default=None,
        help="Static Manhattan plot (non-interactive snapshot)",
    )
    parser.add_argument(
        "-i",
        "--interactive",
        dest="interactive_file",
        metavar="FILE",
        default=None,
        help="Interactive Manhattan plot",
    )
    parser.add_argument("--sep", default="\t", help="Input delimiter (default: tab)")
    parser.add_argument(
        "--chr",
        dest="chr_col",
        default=None,
        metavar="NAME",
        help="Chromosome column (default: auto-detect from bundled formatbook aliases)",
    )
    parser.add_argument(
        "--pos",
        dest="pos_col",
        default=None,
        metavar="NAME",
        help="Position column (default: auto-detect)",
    )
    parser.add_argument(
        "--p",
        dest="p_col",
        default=None,
        metavar="NAME",
        help="P-value column (default: auto-detect; do not use together with --mlog10p)",
    )
    parser.add_argument(
        "--mlog10p",
        dest="mlog10p_col",
        default=None,
        metavar="NAME",
        help="Per-variant -log10(P) column (default: auto-detect if P is absent)",
    )
    parser.add_argument("--skip", type=float, default=5.0, help="Skip variants below -log10(P)")
    parser.add_argument("--width", type=int, default=100)
    parser.add_argument("--height", type=int, default=28)
    parser.add_argument("--ascii", action="store_true", help="Use ASCII fallback")
    parser.add_argument("--sig-level", type=float, default=5e-8)
    parser.add_argument("--ymax", type=float, default=None)
    parser.add_argument(
        "--gtf",
        default=default_gtf_path(),
        help="Protein-coding GRCh37 gene annotation GTF(.gz) for interactive gene track",
    )
    parser.add_argument(
        "--gtf38",
        default=None,
        help="Protein-coding GRCh38 gene annotation GTF(.gz) for interactive gene track",
    )
    return parser


def _resolve_input_and_mode(args: argparse.Namespace, parser: argparse.ArgumentParser) -> tuple[str, str]:
    """Return (path, mode) with mode in {'static', 'interactive'}."""
    pos = args.sumstats
    static = args.static_file
    inter = args.interactive_file
    n = sum(x is not None and str(x).strip() != "" for x in (pos, static, inter))
    if n == 0:
        parser.error("Provide summary statistics: FILE, -s FILE, or -i FILE")
    if n > 1:
        parser.error("Use only one of: positional FILE, -s FILE, or -i FILE")
    if inter is not None:
        return str(inter), "interactive"
    if static is not None:
        return str(static), "static"
    assert pos is not None
    return str(pos), "interactive"


def _app_version() -> str:
    try:
        return version("gwaspeek")
    except PackageNotFoundError:
        return "dev"


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    path, mode = _resolve_input_and_mode(args, parser)

    if mode == "interactive":
        print(f"[gwaspeek v{_app_version()}] Loading {path} ...", file=sys.stderr, flush=True)

    df = load_sumstats(
        path=path,
        sep=args.sep,
        chrom_col=args.chr_col,
        pos_col=args.pos_col,
        p_col=args.p_col,
        mlog10p_col=args.mlog10p_col,
    )
    clean = preprocess_sumstats(df, skip=args.skip)
    use_unicode = not args.ascii

    if mode == "static":
        out = render_manhattan(
            clean,
            width=args.width,
            height=args.height,
            sig_level=args.sig_level,
            ymax=args.ymax,
            unicode=use_unicode,
            y_min=float(args.skip),
        )
        print(out)
        return

    run_interactive_manhattan(
        clean,
        width=args.width,
        height=args.height,
        sig_level=args.sig_level,
        ymax=args.ymax,
        unicode=use_unicode,
        y_min=float(args.skip),
        gtf_path=args.gtf,
        gtf38_path=args.gtf38,
    )


if __name__ == "__main__":
    main()
