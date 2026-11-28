#!/usr/bin/env bash
# Context flags from saved, normalized reference-loss evidence; no caller replay.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 cat <<'HELP'
Usage: dev/reportncrnalosscontext.sh in=cases.tsv out=fresh_prefix t=1 -Xmx1g
  referencefa=scored.fa.gz orientation=sense expectedrecords=N
  contiglengths=optional.tsv endpercent=25

Case TSV columns: family, variant, before, after, comparison, seeds, expected_losses.
Before/after are existing normalized export prefixes; comparison is the saved
comparator prefix. Seeds is the unchanged forward-K17 FASTA (comma-separated
paths combine members of one biological family). The existing reference FASTA
has escaped GFF ID/input_id attributes and gene-sense sequences; no extraction.
Optional contig lengths have columns input_id, seqid, length. All TSVs use tabs.
Missing context stays NA. Classes overlap; tandem-like geometry is not a causal
failure assignment. Outputs are .loci.tsv, .summary.tsv and .definitions.tsv.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.ReportNcrnaLossContext "$@"
