#!/usr/bin/env bash
# Full signed endpoint-error direction report from existing matched-error bins.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Report complete signed endpoint-error masses from an existing member histogram.
Usage: ncrnasignedendpointerrors.sh in=member.histogram.tsv out=fresh.bands.tsv family=RF00174

Input is NcrnaEndpointSummary's family/end/measure/offset_bp/loci table for ONE
exact member. Output separates negative/zero/positive counts within and beyond
0,6,10,16,20 bases. Tails remain intact; no conditioning, clamping or resampling.
Both ends must have the same matched-locus denominator. This reads saved small
tables only; it performs no genome calling, grading, extraction or training.
Use -Xmx to set heap. Limits, sampling, overwrite and append are unsupported.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.NcrnaSignedEndpointErrors "$@"
