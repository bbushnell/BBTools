#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Partition complete census genomes and annotations without changing sequences.'
  echo 'records=TSV loci=TSV genomes=FASTA candidates=FASTA outprefix=NEW_PREFIX partitions=32 window=4139'
  echo 'expectedrecords=N expectedloci=N expectedsearchbases=N expectedcandidatebases=N bind the FULL input census.'
  echo 'Genomic annotation-first subset: genomiconly=t genomicmaxusable=1; candidate inputs are not emitted.'
  echo 'metadataonly=t counts and records eligibility without reading FASTAs; omit genomes= and candidates=.'
  exit 0
fi
exec java -ea -Xmx2g -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=4 --add-modules jdk.incubator.vector -cp "$root/current" prok.NcrnaMitoPartitions "$@"
