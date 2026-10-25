#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Evaluate every member of a complete mitochondrial cluster arm; no filtering or merge cutoff.'
  echo 'in=ARM.fa.gz assignments=NATIVE.tsv consensus=CLUSTERS.fa.gz hbm=CLUSTERS.hbm.gz audit=CONSTRUCTION.tsv.gz'
  echo 'tree=PINNED_TAXTREE lineagerank=phylum arm=LABEL expected=N verified=SHA80_PATHS outprefix=NEW_PREFIX'
  echo 'ref16=FILE ref16id=universal ref18=FILE ref18id=euk18S_universal ref23=FILE ref23id=universal reflsu=FILE reflsuid=lsu_universal_final'
  echo 'lengthbins=200,500,1000,1500,2000,3000,5000 t=16 -Xmx8g'
  echo 'Outputs: .members.tsv.gz, .pairs.tsv.gz, .summary.tsv, .composition.tsv. Unassigned members remain explicit with NA identities.'
  echo 'Unique means give each unique source sequence equal weight after averaging its assigned-copy identities; partial and pseudo strata can overlap.'
  echo 'Summary scopes are ARM, CLUSTER and LINEAGE, with assigned/unassigned denominators; model n= must equal its kept-assignment count.'
  echo 'Pairs include every unordered cluster pair plus diagonals and all four universals. Shorter sequence is query; equal lengths use the left input.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=8g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -cp "$root/current" prok.NcrnaMitoClusterEvaluation "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
