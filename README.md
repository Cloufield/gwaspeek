# gwast

Terminal Manhattan and QQ plots for GWAS summary statistics.

## Install

```bash
pip install -e .
```

## Input requirements

- Tabular file (TSV/CSV) with GWAS summary stats.
- Required columns by default: `CHR`, `POS`, `P`.
- You can override column names with `--chrom-col`, `--pos-col`, `--p-col`.
- `CHR` accepts numeric chromosomes plus `X`, `Y`, `MT` (or `chr`-prefixed forms).

## CLI usage

### Manhattan plot

```bash
gwast manhattan --input sumstats.tsv
```

Common options:

```bash
gwast manhattan \
  --input sumstats.tsv \
  --sep "\t" \
  --width 120 \
  --height 32 \
  --skip 2 \
  --sig-level 5e-8 \
  --ascii
```

### QQ plot

```bash
gwast qq --input sumstats.tsv
```

Example:

```bash
gwast qq --input sumstats.tsv --width 100 --height 28 --ascii
```

## Notes and limitations

- Output is text-mode for terminal readability and portability.
- Plot quality depends on terminal size (`--width`, `--height`).
- Unicode mode gives better visuals; use `--ascii` for limited terminals.
- Very large GWAS files should use `--skip` to reduce overplotting/noise.