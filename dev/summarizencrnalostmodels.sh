#!/usr/bin/env bash
# Reporting only: group previously classified losses by native original model ID.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/summarizencrnalostmodels.sh in=full.loci.tsv,old.loci.tsv panels=full,old out=fresh_prefix
Inputs are completed CompareCobalaminBenchmarkGrades locus tables with measured
same-generation subset membership. Only graded losses without any same-input,
contig and strand clan overlap contribute model names. Gains and overlap-present
losses contribute nothing. Already-retained winners are reported separately.
Writes panel/member counts, original winner rows, a unique removed-name list and
individual already-retained-winner loci. UNION locus counts count evidence rows;
UNION model counts deduplicate names across panels. Does not export a library,
rematch references, classify overlaps again or run a caller.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.SummarizeNcrnaLostModels "$@"
