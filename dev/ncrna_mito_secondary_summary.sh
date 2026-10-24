#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Summarize frozen secondary-member identities without alignments.'
  echo 'in=MEMBERS.tsv.gz out=NEW_SUMMARY.tsv arm=NAME expected=COPIES verified=PINNED_INPUTS.sha80'
  exit 0
fi
exec java -ea -Xmx2g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=2 --add-modules jdk.incubator.vector -cp "$root/current" prok.NcrnaMitoSecondarySummary "$@"
