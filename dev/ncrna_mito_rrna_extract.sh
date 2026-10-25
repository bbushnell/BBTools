#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Extract every annotated mitochondrial rRNA from original GBFF; no SSU classification.'
  echo 'in=GBFF roster=TSV fna=FULL_GENOMIC_FASTA expectedrecords=N outprefix=PATH t=16 zl=9'
  echo 'Outputs: .candidates.fa.gz .genomes.fa.gz .loci.tsv.gz .records.tsv .summary.tsv'
  echo 'Retains partial/pseudo features and all roster genomes. Unsupported geometry stays in the ledger.'
  echo 'fna= supplies accession-bound full sequences for CONTIG-only records; original sequence provenance is recorded.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=8g --mode=fixed "$@";then
  echo 'Invalid Java launcher arguments' >&2
  exit 1
fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
command=(java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=16 -cp "$root/current" prok.NcrnaMitoRrnaExtract "$@")
printf '%q ' "${command[@]}" >&2
printf '\n' >&2
exec "${command[@]}"
