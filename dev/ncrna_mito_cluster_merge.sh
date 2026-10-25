#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Report single-linkage groups at every specified identity cut from saved evaluation outputs; no alignments or cutoff selection.'
  echo 'in=EVALUATION.members.tsv.gz pairs=EVALUATION.pairs.tsv.gz tree=PINNED_TAXTREE verified=SHA80_PATHS arm=NAME expected=N expectedmodels=K out=NEW_TABLE outtree=NEW_TREE'
  echo 'cuts=0.60,0.65,0.70,0.75,0.80,0.85,0.90 t=1 -Xmx4g'
  echo 'One row per top-three group at each cut; fewer rows when fewer than three groups exist. Bridges are transitive single linkage.'
  echo 'Rank by member copies, unique sequences, then lexical group key. Dominant product/phylum/class/best-CM labels use copy counts; all modal ties and support counts are reported.'
  echo 'Unassigned members remain in arm/unassigned denominators. Labels describe source annotation and host taxonomy, not inferred RNA ancestry.'
  echo 'The merge tree records actual maximum-spanning union edges, descending identity then lexical endpoint names; it reproduces all single-linkage cuts.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=4g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 -cp "$root/current" prok.NcrnaMitoClusterMerge "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
