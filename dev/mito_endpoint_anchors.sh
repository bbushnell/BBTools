#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Capture shipped raw mitochondrial endpoints: in=FLANKS.fa.gz ledger=PILOT.tsv.gz expected=N outprefix=NEW_PREFIX t=4'
  exit 0
fi
exec java -ea -Xmx8g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.MitoEndpointAnchors "$@"
