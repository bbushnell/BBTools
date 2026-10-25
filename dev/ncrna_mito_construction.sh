#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Audit accepted mitochondrial evidence and export locked m16S/m18S annotation cohorts.'
  echo 'in=SOURCE_CANDIDATES loci=EVIDENCE_LOCI hits=EVIDENCE_HITS associations=EVIDENCE_ASSOCIATIONS pairs=EVIDENCE_UNIVERSAL_PAIRS genomes=EVIDENCE_GENOMES verified=SHA80_PATHS outprefix=NEW_PREFIX'
  echo 'expectedloci=37716 expectedgenomes=18236 expectedhits=89876 expectedgroups=89773 expected16=4831 expected18=6586 expectedunsupported=8672 t=1 zl=9 -Xmx2g'
  echo 'Partial/pseudo/discordant members are retained. CM-only evidence is audit-only. No searches or production changes.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=2g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 -cp "$root/current" prok.NcrnaMitoConstruction "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
