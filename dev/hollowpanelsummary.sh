#!/bin/bash
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
if [[ $# == 0 || "$1" == --help || "$1" == -h ]]; then
  echo "Summarize a complete matched solid/Hollow coordinate-audit panel."
  echo "in=samples.tsv root=audit_root out=paired.tsv details=pooled.tsv genomes=100"
  echo "Requires seven settings, two geometries and caps20/100 for every genome."
  echo "pairmanifest=t instead accepts explicit baseline/treatment pairs with one common cap."
  echo "Its in= TSV header is: genome cap baseline_arm baseline_dir treatment_arm treatment_dir (tabs)."
  echo "Each directory contains summary.tsv and, for cases=, suffixes.tsv; root= is unused."
  echo "Arm labels must match the original audit labels. One row per genome; defaults: pairmanifest=f, genomes=100."
  echo "Preserves per-genome regressions. Nonexact suffixes are not inferred MSDI counts."
  echo "cases=file optionally joins all per-ID adverse reports, including Hollow-only failures."
  echo "Missing opposite rows mean exact OR unextended; their counts are reported as -1, not zero."
  echo "Circular-compatible denotes sequence agreement, not verified circular topology."
  echo "test=t runs deterministic counter checks. Outputs must be new."
  exit
fi
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowPanelSummary "$@"
