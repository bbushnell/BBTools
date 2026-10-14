#!/bin/bash
# Diagnostic oracle or NN rescue of saved same-K raw unitigs; not production integration.
# in=raw.fa candidates=candidates.tsv labels=labels.tsv k=64 out=merged.fa
# outbaseline=baseline.fa outpaths=paths.tsv mincontig=500
# OR replace labels= with features=features.tsv net=model.bbnet cutoff=0.8 outscore=scores.tsv
# NN mode scores ALL candidate rows, never reads truth labels, and accepts score>=cutoff.
# It uses scalar inference to match validation and is limited to the K24–127 training domain.
# Accepted decisions with unique reciprocal endpoint choices may form chains.
# Closed cycles are left unmerged. Filtering occurs after merging, identically
# for both FASTA outputs. No read tables or references are loaded by this step.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKJoinMerge "$@"
