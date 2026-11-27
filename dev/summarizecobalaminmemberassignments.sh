#!/usr/bin/env bash
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/summarizecobalaminmemberassignments.sh in=baseline.acceptance_calls.tsv \
 in2=lower.acceptance_calls.tsv grades=baseline,lower,actual_global \
 reference=references.tsv raw=baseline.overlap.raw.tsv mapping=members.tsv \
 evidence=sources.tsv choice=RF01482.band_choice.tsv \
 turnover=RF01482.member_turnover.tsv losses=global_loss.details.tsv \
 member=RF01482 nativepass=.66 expectedinputs=10000 expectedreferences=21524 \
 expectedmasked=2323 expectedmemberlost=151 expectednosame=126 \
 out=fresh_prefix t=1 -Xmx2g -ea -eoom

Cross-tabs frozen emitted members against two independent annotation axes:
stored representative and actual protected-prior/confident-witness support.
Multiple supported members remain AMBIGUOUS; valid evidence without confident
support remains UNKNOWN. Missing or malformed required evidence fails visibly.
Unmatched calls have NO_MATCHED_REFERENCE; no overlap-derived label is invented.

Reapplies one completed descriptive member-band choice to saved scores and
reconciles every split against its turnover table. Separately reports the
actual global-filter losses with literal NO_SAME_CM_MEASUREMENT status.
The all-match control cohorts include gained references. Masked calls are
excluded and accounted for. Cohorts overlap and must not be added together.

Exact original input-qualified reference IDs and raw query assignments bind
the metadata. GFFs and explicit witness/event sidecars are read; FASTA contents
are not read. No caller, CM rescoring, rematching, fit or policy adoption.

Writes .matrix.tsv and .cohorts.tsv for delivery; .details.tsv and
.annotations.tsv retain native audit evidence including full support sets.
Annotation agreement is descriptive and is not biological subtype accuracy.
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
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError \
 -cp "${DIR}current/" prok.SummarizeCobalaminMemberAssignments "$@"
