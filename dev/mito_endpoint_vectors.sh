#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'in=FLANKS.fa.gz anchors=ANCHORS.tsv.gz ledger=SOURCE_WINDOW_LEDGER.tsv.gz consensus=CONSENSUS.fa consensussha80=PIN expected=256 outprefix=NEW_PREFIX startmin=N startmax=N stopmin=N stopmax=N startinside=N stopinside=N t=1'
  echo 'Explicit k9 placement and measured correction ranges are required. TRAIN-only tables; source-host splits are preserved. No network fit.'
  echo 'Optional net=MODEL.bbnet tables=SAVED.tables.tsv.gz evaluates the actual sequential scorer on validation anchors instead of writing vectors. This is a development harness, not caller accuracy.'
  exit 0
fi
exec java -ea -Xmx8g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.MitoEndpointVectors "$@"
