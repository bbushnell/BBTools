#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Export an exact single-linkage component with full copy/flag conservation; no taxon filtering.'
  echo 'in=SOURCE.fa.gz models=MODELS.fa.gz members=MEMBERS.tsv.gz pairs=PAIRS.tsv.gz arm=NAME target=CANONICAL_KEY group=NAME cut=.70 expected=N expectedmodels=N selected=N outprefix=NEW_PREFIX verified=INPUTS.sha80'
  exit 0
fi
exec java -ea -Xmx2g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=2 --add-modules jdk.incubator.vector -cp "$root/current" prok.NcrnaMitoGroupExport "$@"
