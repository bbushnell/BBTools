#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'in=FULL_FLANKS.fa.gz ledger=FULL_SOURCE.tsv.gz assignments=FULL_ASSIGNMENTS.tsv.gz consensus=SHIPPED.fa consensussha80=PIN outprefix=FRESH_PREFIX startinside=9 stopinside=9 t=1'
  echo 'Keeps16110members/15084usable anchors, builds TRAIN-only k9 tables and both radius128/256 vector sets in one pass. Reports zero-positive blocks explicitly. No network fit.'
  echo 'Evaluation only: add net=EXPORTED.bbnet tables=TRAIN_ONLY.tsv.gz radius=192 contextgenomes=ORIGINAL_GENOMES.fa.gz. Expands artificial crops with verified source context while retaining raw anchors. Scores all3019usable validation anchors; does not generate vectors or install a caller network.'
  echo 'Optional thresholds=-Infinity,0,0.1,0.5,Infinity evaluates shared absolute-score gates with cached sequential geometry. Include both infinity controls; scores are not probabilities.'
  exit 0
fi
exec java -ea -Xmx8g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 --add-modules jdk.incubator.vector -cp "$root/current" test.prok.MitoEndpointTrainingVectors "$@"
