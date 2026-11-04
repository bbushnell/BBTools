#!/usr/bin/env bash
# Reconstruct the full declared cohort; emit a bounded cost sample, never train.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'ssu=SSU.fa.gz insect=INSECT.fa.gz metadata=SOURCE.tsv.gz genomes=GENOMIC.fa.gz expectedssu=N expectedinsect=N sample=N flank=N maxspan=N seed=N outprefix=NEW_PREFIX t=4'
  echo 'Optional exact-roster recovery: enrollment=COVERAGE_MEMBERS.tsv.gz insectinput=7050, with insect=the complete nohit FASTA. Only declared CONSTRUCTION IDs are selected; PANEL_INCLUDED rows are excluded.'
  exit 0
fi
exec java -ea -Xmx8g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.MitoEndpointCorpus "$@"
