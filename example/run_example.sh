#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
# Default to the real gwaslab sample dataset requested by user.
# You can override via: INPUT=/path/to/file bash example/run_example.sh
INPUT="${INPUT:-/home/yunye/work/gwaslab/examples/0_sample_data/t2d_bbj.txt.gz}"

echo "== Manhattan (ASCII) =="
PYTHONPATH="${ROOT_DIR}/src" python -m gwast.cli manhattan \
  --input "${INPUT}" \
  --chrom-col CHR \
  --pos-col POS \
  --p-col P \
  --skip 3 \
  --width 80 \
  --height 22 \
  --ascii

echo
echo "== QQ (ASCII) =="
PYTHONPATH="${ROOT_DIR}/src" python -m gwast.cli qq \
  --input "${INPUT}" \
  --chrom-col CHR \
  --pos-col POS \
  --p-col P \
  --skip 3 \
  --width 80 \
  --height 22 \
  --ascii
