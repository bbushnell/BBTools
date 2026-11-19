#!/usr/bin/env bash
# Projects already drawn held-out trials; no alignment, fitting or table teaching.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Project saved 69-input ncRNA endpoint trials into the 24-input inside-K9 layout.
Usage: projectncrnapositionalinputs.sh in=evaluation.tsv rows=evaluation.rows.tsv \
 contexts=original.fa roster=original.genomes.tsv ref=selected.consensus.fa \
 table=retained.5prime.table.tsv member=RF00174 end=5prime expected=13404 \
 out=fresh24.tsv

Requires radius10 and the retained fully-inside K9 table for the declared end.
Reuses the original sidecar directly. Validates source, assembly, true spans,
selected model names/indices and full trial conservation. Copies the common
length/GC/identity features and one-hot labels verbatim; recomputes only 21 scores.
The original context FASTA is read once, retaining only needed locus bases.
No trial resampling, new model selection, alignment, table teaching or training.
Outputs are always fresh. A partial file from a failed process is not usable.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.ProjectNcrnaPositionalInputs "$@"
