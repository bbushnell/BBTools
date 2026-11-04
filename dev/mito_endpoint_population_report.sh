#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'evaluation=EVALUATION.tsv.gz enrollment=ENROLLMENT.tsv.gz tree=FROZEN_TAXONOMY.tsv.gz rules=LINEAGES.tsv outprefix=FRESH_PREFIX [gates=GATES.tsv.gz]'
  echo 'Reports exact and both/per-end within0/1/2/5/10/20, full signed errors and lineage conservation. Optional gates report aggregate nonregression and individual gains/losses.'
  exit 0
fi
exec java -ea -Xmx8g -XX:ActiveProcessorCount=2 -XX:+ExitOnOutOfMemoryError --add-modules jdk.incubator.vector -cp "$root/current" test.prok.MitoEndpointPopulationReport "$@"
