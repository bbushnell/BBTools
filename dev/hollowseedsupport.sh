#!/bin/bash
# Read-only evidence diagnostic; never edits seeds, counts or extensions.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
if [[ $# == 0 || "$1" == --help || "$1" == -h ]]; then
  echo "Measure original-read solid support for saved single-base seed changes."
  echo "in=original.fq.gz seeds=original.fa queries=extended.fa out=changes.tsv stats=summary.tsv"
  echo "k=96 (64 or96), support=2, t=N; outputs must be new."
  echo "Original seeds and proposed outputs must contain identical ordered IDs."
  echo "At most one substitution per retained seed is allowed; unchanged seeds are skipped."
  echo "Counts all original reads without trimming, merging, sampling or prefiltering."
  echo "Queries every solid window containing the change, wholly inside the original seed."
  echo "Reports depth extrema and supported-window counts for original and proposed bases."
  echo "No reference truth, filtering decision or sequence output. Appended bases are excluded."
  echo "test=t in=existing_small.fa runs deterministic in-memory geometry/count tests."
  exit
fi
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowSeedSupport "$@"
