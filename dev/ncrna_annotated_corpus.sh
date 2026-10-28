#!/usr/bin/env bash
# Extract an Rfam population using the survey's source-locus definition.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 ]];then
  printf '%s\n' 'Usage: ncrna_annotated_corpus.sh roster=survey.genomes.tsv family=RF00174 expectedgenomes=39979 outprefix=corpus t=8 [limit=300 flank=200 expectedloci=91499]'
  exit 0
fi
exec java -ea --add-modules jdk.incubator.vector -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=8 -Xmx8g -cp "$root/current" prok.NcrnaAnnotatedCorpus "$@"
