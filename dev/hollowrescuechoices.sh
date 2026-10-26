#!/bin/bash
# Read-only post-hoc evidence; never changes extension decisions or sequences.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
if [[ $# == 0 || "$1" == --help || "$1" == -h ]]; then
  echo "Reconstruct all16 Hollow ingress/next-base counts for saved right-extension seed edits."
  echo "in=original.fq.gz seeds=original.fa queries=extended.fa out=choices.tsv stats=summary.tsv"
  echo "k=192 support=4 t=N; outputs must be new."
  echo "Counts original reads AND original seeds exactly once, without filtering."
  echo "Queries never contribute evidence; no reference truth is consulted."
  echo "Ordered unpaired FASTA comparisons must retain seed length and at most one seed substitution."
  echo "Checkpoint=edited_position+k-32; only retained right-extension rescue geometry is accepted."
  echo "Counts use query orientation: count_AC means left ingressA and next baseC."
  echo "supported_competitor excludes all four successors of the proposed ingress."
  echo "chosen_not_maximum flags reconstruction/input disagreement; this is not an execution trace."
  echo "test=t in=existing_small.fa runs deterministic in-memory geometry/count tests."
  exit
fi
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=32g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowRescueChoices "$@"
