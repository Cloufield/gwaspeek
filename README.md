# gwaspeek

Terminal **interactive** Manhattan plots for GWAS summary statistics: pan and zoom in a real TTY, inspect regions, and optionally show a protein-coding **gene track** when the view is a single chromosome and spans **≤ 1 Mb**.

**Also:** one-shot **static** renders (`-s`) for logs or CI; column names **auto-detected** from bundled [formatbook](https://github.com/Cloufield/formatbook)-style aliases; **Unicode** or **`--ascii`**.

## Installation

```bash
pip install -e .
gwaspeek --help
```

## Try it

After install, from a directory that contains the file:

```bash
gwaspeek t2d_bbj_p1e-5.txt.gz
```

Interactive mode is the default whenever you pass exactly one path as **`FILE`** or **`-i FILE`**. Use a real TTY when you can. Keys: **`A`** / **`D`** pan, **`W`** / **`S`** zoom out / in, **`h`** help, **`q`** quit — full list under [Interactive keys](#interactive-keys).

### Static preview (same data)

Text snapshot (non-default):

```bash
gwaspeek -s t2d_bbj_p1e-5.txt.gz --width 120 --ascii
```

Sample output:

```text
Manhattan 1:17276978-X:153269341                                                                                gwaspeek
 205.0│
      │
      │
      │
      │                                                                     ●
      │                                                                     ●
 155.0│                                                                     ●
      │                                                                     ●
      │
      │
      │                                                                     ●
      │
 105.0│                                                                     ●
      │                                                                     ●
      │                                                           ●         ●
      │                                        ○                  ●         ●
      │                                        ○          ●       ●         ●
      │                                        ○          ●       ●         ●
   55.0│                                        ○          ●                                                           ○
      │                                        ○          ●       ●         ○                                         ○
      │                         ●              ○                  ●        ○●                                         ○
      │                   ●     ○              ○          ●     ○ ●        ○○                         ●
      │         ○         ●●    ●     ○   ●    ○○   ○   ● ●  ○  ○ ●  ● ●  ○○○     ○○   ●● ●    ●●  ○ ●●   ○ ●        ○○
   5.0│╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌○
      +───+─────────+───────+───────+──────+──────+─────+─────+─────+─────+────+────+───+───+──+───+──+──+─+─+─+───+────
          1         2       3       4      5      6     7     8     9    10   11   12  13  14 15  16 17 18192022   X
```

## How to invoke

Give **exactly one** input path:

| Form | Mode |
|------|------|
| `gwaspeek FILE` | Interactive (same as `-i`) |
| `gwaspeek -i FILE` | Interactive |
| `gwaspeek -s FILE` | Static |

More than one path (e.g. positional plus `-s`) is an error.

Common flags (both modes unless noted):

```bash
gwaspeek sumstats.tsv \
  --sep "\t" \
  --width 120 \
  --height 32 \
  --skip 2 \
  --sig-level 5e-8
```

Add **`-s`** before the path for a static render instead of interactive.

## Input files

Tabular GWAS summary stats (TSV/CSV or other delimiter via **`--sep`**).

**Per row you need:** chromosome, base-pair position, and either **`P`** or **`-log10(P)`** (names are auto-detected; override with **`--chr`**, **`--pos`**, **`--p`**, **`--mlog10p`**). Do not pass both **`--p`** and **`--mlog10p`**. Internally rows are normalized to `CHR`, `POS`, `P`, and derived `mlog10p`.

**Chromosome tokens** understood by the parser include numeric `1`–`22`, `X`, `Y`, `MT`, and `chr`-prefixed forms (e.g. `chr1`, `chrX`, `chrM`).

## CLI reference

| Option | Type | Default | Description |
|---|---|---|---|
| `FILE` (positional) | path | — | Same as **`-i FILE`** (interactive). |
| `-i`, `--interactive FILE` | path | — | Interactive viewer. |
| `-s`, `--static FILE` | path | — | One-shot static plot. |
| `--sep STR` | string | tab (`\t`) | Delimiter. |
| `--chr NAME` | string | auto | Chromosome column. |
| `--pos NAME` | string | auto | Position column. |
| `--p NAME` | string | auto | P-value column (not with **`--mlog10p`**). |
| `--mlog10p NAME` | string | auto | `-log10(P)` column (not with **`--p`**). |
| `--skip FLOAT` | float | `5.0` | Drop variants below this `-log10(P)`; also y-axis floor. |
| `--width INT` | int | `100` | Width (chars): static size; interactive initial/fallback. |
| `--height INT` | int | `28` | Height (lines): static size; interactive initial/fallback. |
| `--ascii` | flag | off | ASCII drawing instead of Unicode. |
| `--sig-level FLOAT` | float | `5e-8` | Genome-wide line in `P` space. |
| `--ymax FLOAT` | float | auto | Y-axis max in `-log10(P)`. |
| `--gtf PATH` | path | bundled GRCh37 | GTF (`.gz` ok); interactive gene track. |
| `--gtf38 PATH` | path | — | GRCh38 GTF (`.gz` ok); interactive mode uses bundled gene-only GTF when omitted. |
| `-h`, `--help` | flag | off | Help. |

## Interactive keys

With **`gwaspeek FILE`** or **`gwaspeek -i FILE`**:

| Key | Action |
|-----|--------|
| `A` / `D` | Pan |
| `W` / `S` | Zoom out / in |
| `l` | Toggle lead-variant labels |
| `t` | Cycle gene track: `37` → `38` → `off` |
| `v` | Toggle variants-in-view list |
| `h` | Toggle help |
| `r` | Reset to full genome |
| `q` | Quit |

Gene track appears only on a **single-chromosome** view with span **≤ 1 Mb**. **`--gtf`** / **`--gtf38`** apply to interactive mode only; static **`-s`** output ignores them.

## Examples

```bash
# Interactive (default): positional or -i
gwaspeek tests/fixtures/sumstats_small.tsv
gwaspeek -i tests/fixtures/sumstats_small.tsv

# Static
gwaspeek -s tests/fixtures/sumstats_small.tsv

# Explicit columns (static shown; same flags work for interactive)
gwaspeek -s sumstats.tsv --chr CHROM --pos BP --p PVAL

# Precomputed -log10(P)
gwaspeek -s sumstats_mlog10p.tsv --chr CHR --pos POS --mlog10p LOG10P

# Narrow terminals
gwaspeek -s sumstats.tsv --ascii
```

## Tips

- Prefer a **real TTY** for interactive mode.
- Dense sumstats: raise **`--skip`** to thin points and reduce clutter.
- Unicode is usually clearer; **`--ascii`** when the terminal cannot draw box drawing / braille well.

## For contributors

```bash
pip install -e ".[dev]"
pytest
```

Pipeline: **`gwaspeek.cli`** parses args and loads data with **`gwaspeek.io`**, normalizes with **`gwaspeek.preprocess`**, then **`gwaspeek.manhattan.render_manhattan`** (static) or **`gwaspeek.interactive.run_interactive_manhattan`** (interactive). Terminal drawing lives in **`gwaspeek.terminal_canvas`**. Tests live under **`tests/`** (see `tests/conftest.py` for `PYTHONPATH` / `src` layout).

## License

MIT.

## Repository

[https://github.com/Cloufield/gwast](https://github.com/Cloufield/gwast)
