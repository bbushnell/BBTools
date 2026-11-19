#!/usr/bin/env bash
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/ncrnaboundaryinvocationsummary.sh in=invocations.tsv out=fresh.tsv expectedinputs=1000 t=1 -Xmx256m -eoom
Summarizes invocation-only boundary traces. Candidate custom-refiner calls,
returned events and retained emitted calls are separate counts. These are not
matched-reference loci. Unknown provenance remains separate from known holds.
Validates the distinct header, source declarations and candidate/contig totals.
No calling, grading, inference or score replay occurs.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.NcrnaBoundaryInvocationSummary "$@"
