#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'in=FLANKS.fa.gz ledger=SOURCE_WINDOWS.tsv.gz ssu=SSU.fa.gz insect=NOHIT.fa.gz consensus=SHIPPED.fa pilotanchors=ORIGINAL256.tsv.gz expected=256|16110 outprefix=NEW_PREFIX t=1 passes=1'
  echo 'passes>1 is a single-JVM timing repetition only for expected=256,t=1; checks exact repeated results, writes each member once and reports actual repeated alignment work.'
  echo 'Shared k9 top8, maximum shorter-query identity, one genomic fit per assigned model. Every member retained. No seed gate, all-pairs matrix, training or production change.'
  exit 0
fi
exec java -ea -Xmx8g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.MitoEndpointAssignments "$@"
