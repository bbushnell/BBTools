#!/bin/bash
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
if [[ $# == 0 || "$1" == --help || "$1" == -h ]]; then
  echo "Find literal sequence matches to the outer32bp anchors of saved long-kmer queries."
  echo "in=depth.tsv ref=genome.fa ends=251,257 out=hits.tsv details=summary.tsv label=query"
  echo "Searches both strands at the saved span; no hashing, alignment or circular wrapping."
  echo "Coordinates start0/end0 are zero-based inclusive. Next bases are query-oriented."
  echo "Reports interior differences and full-word equality; both strands count as oriented hits."
  echo "successors=f  With t, also scan all four one-base-shifted A/C/G/T successor keys."
  echo "Appends source_end1, query_role and candidate_base columns only with successors=t."
  echo "Candidate keys are derived from saved words, never from reference truth."
  echo "ref may be unpaired FASTQ (interleaved=f); differences are mismatches, not edit events."
  echo "paths=paths.fa flanks=flanks.tsv  Optionally compare covered read flanks to saved paths."
  echo "Requires exact query IDs and word/path agreement. Clips to both read and path bounds."
  echo "queryfasta=t scaffold=ID k=224: in is a saved FASTA assembly, not a depth table."
  echo "Extracts windows ending at ends= directly from that scaffold; no count table is built."
  echo "scaffold matches the first header token; full headers remain the paths= identity."
  echo "Flank differences are coordinate mismatches, not alignments or error-event counts."
  echo "Input-record ordinals distinguish duplicate names; standard hit/summary output is unchanged."
  echo "test=t runs deterministic matching checks. Outputs must be new."
  exit
fi
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowAnchorOrigins "$@"
