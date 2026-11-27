#!/usr/bin/env bash
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/summarizecobalaminmemberendpoints.sh in=accepted.details.tsv \
 in2=accepted.annotations.tsv witnesses=lower.normalized.witnesses.tsv \
 metrics=lower.normalized.metrics.tsv capture=SAVED_LOWER_CAPTURE \
 expectedkept=1414 expectedselected=1126 out=fresh.tsv t=1 -Xmx1g -ea -eoom

Reports both-ends-exact counts and fractions for the selected RF01482 calls
on RF00174-only supported references, and for RF00174 calls from that SAME
saved capture. The main-only comparator and the broader any-main-supported
comparator include all matching calls in that capture, without a new filter.
The latter includes the former; do not add these comparator populations.

Reads closed assignment/support IDs and saved normalized d5/d3 values. Both
ends exact means d5=0 AND d3=0. Reconciles the complete capture's matched and
both-exact counts and every selected member/model/reference binding. No caller,
rematching, support regeneration, endpoint refinement, fit or matrix rerun.
HELP
 exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ];do
 cd "$(dirname "$SCRIPT")"
 nextScript="$(readlink "$(basename "$SCRIPT")")"
 SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError \
 -cp "${DIR}current/" prok.SummarizeCobalaminMemberEndpoints "$@"
