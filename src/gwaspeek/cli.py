from __future__ import annotations

import argparse
import sys
from gwaspeek.interactive import default_gtf_path, run_interactive_manhattan
from gwaspeek.io import load_sumstats
from gwaspeek.manhattan import render_manhattan, stdout_color_supported
from gwaspeek.plot_state import prepare_plot_dataset
from gwaspeek.preprocess import preprocess_sumstats
from gwaspeek.versioning import package_version


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Keep epilog/examples literal; append default values to option lines."""


def build_parser() -> argparse.ArgumentParser:
    ver = package_version()
    parser = argparse.ArgumentParser(
        prog="gwaspeek",
        description=(
            f"gwaspeek {ver} — draw GWAS Manhattan plots in the terminal: static snapshot (-s) or "
            f"full-screen interactive pan/zoom. Show the installed version with -v / --version."
        ),
        formatter_class=_HelpFormatter,
        epilog=(
            "Examples:\n"
            "  gwaspeek tests/fixtures/sumstats_small.tsv\n"
            "  gwaspeek -i tests/fixtures/sumstats_small.tsv\n"
            "  gwaspeek -s tests/fixtures/sumstats_small.tsv\n"
            "\n"
            f"Version: gwaspeek {ver}  (also: gwaspeek -v  or  gwaspeek --version)\n"
        ),
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"gwaspeek {package_version()}",
        help="Print gwaspeek version and exit",
    )
    parser.add_argument(
        "sumstats",
        nargs="?",
        metavar="FILE",
        help="GWAS summary statistics file (TSV/CSV). Same as -i FILE when this is the only input.",
    )
    parser.add_argument(
        "-s",
        "--static",
        dest="static_file",
        metavar="FILE",
        default=None,
        help="Print one Manhattan frame to stdout and exit (non-interactive)",
    )
    parser.add_argument(
        "-i",
        "--interactive",
        dest="interactive_file",
        metavar="FILE",
        default=None,
        help="Open the interactive Manhattan viewer (explicit form of positional FILE)",
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Field delimiter in the input file (default: %(default)r)",
    )
    parser.add_argument(
        "--chr",
        dest="chr_col",
        default=None,
        metavar="NAME",
        help="Chromosome column name (auto-detect from bundled formatbook aliases if omitted)",
    )
    parser.add_argument(
        "--pos",
        dest="pos_col",
        default=None,
        metavar="NAME",
        help="Genomic position column name (auto-detect if omitted)",
    )
    parser.add_argument(
        "--p",
        dest="p_col",
        default=None,
        metavar="NAME",
        help="P-value column name (auto-detect; do not combine with --mlog10p)",
    )
    parser.add_argument(
        "--mlog10p",
        dest="mlog10p_col",
        default=None,
        metavar="NAME",
        help="Per-variant -log10(P) column name (auto-detect if P is absent)",
    )
    parser.add_argument(
        "--skip",
        type=float,
        default=5.0,
        help="Hide variants with -log10(P) below this threshold (also used as the y-axis floor)",
    )
    parser.add_argument("--width", type=int, default=100, help="Terminal width in characters for the plot frame")
    parser.add_argument("--height", type=int, default=28, help="Terminal height in lines for the plot frame")
    parser.add_argument(
        "--ascii",
        action="store_true",
        help="Use ASCII drawing characters instead of Unicode box drawing and glyphs",
    )
    parser.add_argument(
        "--build",
        choices=["37", "38"],
        default="37",
        help="Reference assembly for cytoband lengths and cumulative genome layout",
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable ANSI colors (interactive status/footer and static plot chromosome colors)",
    )
    parser.add_argument(
        "--sig-level",
        type=float,
        default=5e-8,
        metavar="P",
        help="Genome-wide significance P-value; drawn as a horizontal threshold in -log10(P) space",
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        metavar="L",
        help="Clamp the plot's maximum -log10(P) to this value (omit for data-driven ceiling)",
    )
    parser.add_argument(
        "--gtf",
        default=default_gtf_path(),
        help="GRCh37 protein-coding gene GTF(.gz) for the interactive gene track. Default: %(default)s",
    )
    parser.add_argument(
        "--gtf38",
        default=None,
        help="GRCh38 protein-coding gene GTF(.gz) for the interactive gene track when build is 38. Default: %(default)s",
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
    return package_version()


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
    prepared = prepare_plot_dataset(clean, build=args.build)
    use_unicode = not args.ascii

    if mode == "static":
        chr_color = stdout_color_supported() and not args.no_color
        out = render_manhattan(
            clean,
            width=args.width,
            height=args.height,
            sig_level=args.sig_level,
            ymax=args.ymax,
            unicode=use_unicode,
            y_min=float(args.skip),
            prepared=prepared,
            color=chr_color,
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
        build=args.build,
        color=not args.no_color,
    )


if __name__ == "__main__":
    main()
