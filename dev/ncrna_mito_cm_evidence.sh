#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Preserve mitochondrial pilot CM evidence and physical-locus links without admitting truth.'
  echo 'selection=TSV loci=TSV models=TSV runsroot=DIR members=TSV pairs=TSV verified=SHA80_PATHS outprefix=NEW_PREFIX'
  echo 'expectedgenomes=16 expectedloci=22 expectedmodels=9 t=1 zl=9 -Xmx2g'
  echo 'For partitioned scans: partitions=N membership=TSV; runsroot contains MODEL_UNIVERSE_pINDEX. Omit both for unsplit scans.'
  echo 'Restricted genomic discovery also needs genomicroot=DIR genomicpartitions=N genomicmembership=TSV genomiceligibility=TSV expectedgenomicgenomes=N expectedgenomicbases=N.'
  echo 'Candidates retain the full census; excluded genomic evidence is NA/NOT_SEARCHED and reported separately.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=2g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 -cp "$root/current" prok.NcrnaMitoCmEvidence "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
