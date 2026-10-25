#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Select deterministic annotated/no-rRNA mitochondrial CM pilot targets; preserve all selected labels.'
  echo 'records=TSV loci=TSV genomes=FASTA candidates=FASTA outprefix=NEW_PREFIX expectedrecords=N expectedloci=N expectedkinetoplast=N'
  echo 'pergroup=8 maxbases=5m window=N t=1 zl=9 -Xmx2g. window is the maximum W of the pinned CM roster.'
	echo 'selection=quartile (default) uses pergroup quotas; selection=longtail chooses longest annotated/unannotated, longest member of each special stratum, and longest-candidate genome (at most7); pergroup is unused.'
	echo 'population.tsv reports full-census bases and both-strand Mb, for fixed -Z/--FZ across future shards. No search is run by this tool.'
  echo 'Circular targets append min(window-1,length-1) prefix bases; original genomes and coordinate map are retained.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=2g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 -cp "$root/current" prok.NcrnaMitoPilot "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
