#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Refine one mitochondrial consensus from the complete cohort; duplicate copies retain votes.'
  echo 'in=COHORT.fa.gz seeds=DIAGNOSTIC_CLUSTERS.fa name=MODEL expected=COPIES outprefix=NEW_PREFIX'
  echo 'maxrounds=20 finaltrim=0.5 roundtrim=0 t=16 -Xmx8g; no identity exclusion; nonconvergence fails loudly.'
  echo 'roundtrim=0.5 enables the explicit secondary experiment: half-max depth trimming every round.'
  echo 'tolerancestable=t retains exact termination, but emits the last round if length change<=2 and native Quantum identity>=0.999.'
  echo 'Tolerance requires equal nonzero roundtrim/finaltrim; output is labelled TOLERANCE_STABLE, never exact-stable.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=8g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -cp "$root/current" prok.NcrnaMitoConsensus "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
