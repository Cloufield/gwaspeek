from __future__ import annotations

import argparse

from gwast.io import load_sumstats
from gwast.manhattan import render_manhattan
from gwast.preprocess import preprocess_sumstats
from gwast.qq import render_qq


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="gwast", description="Terminal GWAS plotting")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_shared(p: argparse.ArgumentParser) -> None:
        p.add_argument("--input", "-i", required=True, help="Input summary stats file")
        p.add_argument("--sep", default="\t", help="Input delimiter (default: tab)")
        p.add_argument("--chrom-col", default="CHR")
        p.add_argument("--pos-col", default="POS")
        p.add_argument("--p-col", default="P")
        p.add_argument("--skip", type=float, default=2.0, help="Skip variants below -log10(P)")
        p.add_argument("--width", type=int, default=100)
        p.add_argument("--height", type=int, default=28)
        p.add_argument("--ascii", action="store_true", help="Use ASCII fallback")

    p_m = sub.add_parser("manhattan", help="Render Manhattan plot")
    add_shared(p_m)
    p_m.add_argument("--sig-level", type=float, default=5e-8)
    p_m.add_argument("--ymax", type=float, default=None)

    p_q = sub.add_parser("qq", help="Render QQ plot")
    add_shared(p_q)
    p_q.add_argument("--ymax", type=float, default=None)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    df = load_sumstats(
        path=args.input,
        sep=args.sep,
        chrom_col=args.chrom_col,
        pos_col=args.pos_col,
        p_col=args.p_col,
    )
    clean = preprocess_sumstats(df, skip=args.skip)
    use_unicode = not args.ascii

    if args.command == "manhattan":
        out = render_manhattan(
            clean,
            width=args.width,
            height=args.height,
            sig_level=args.sig_level,
            ymax=args.ymax,
            unicode=use_unicode,
        )
    elif args.command == "qq":
        out = render_qq(
            clean,
            width=args.width,
            height=args.height,
            ymax=args.ymax,
            unicode=use_unicode,
        )
    else:
        parser.error(f"Unknown command: {args.command}")
        return
    print(out)


if __name__ == "__main__":
    main()
