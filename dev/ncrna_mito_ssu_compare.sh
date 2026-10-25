#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Compare every mitochondrial rRNA candidate to the universal16S and euk18S_universal consensuses.'
  echo 'in=FASTA ref16=FASTA ref18=FASTA outprefix=PATH expected=N t=16 zl=9 -Xmx8g'
  echo 'Outputs .pairs.tsv.gz .members.tsv.gz .unique.fa.gz .summary.tsv; no biological classification.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=8g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=16 -cp "$root/current" prok.NcrnaMitoSsuCompare "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
