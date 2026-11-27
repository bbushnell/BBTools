#!/usr/bin/env bash
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/selectcobalaminclanacceptance.sh in=clan.scores.tsv
 in2=analysis.event_assignments.tsv raw=caller.work.acceptance.tsv
 baseline=ordinary.normalized capture=lower.normalized sources=sources.tsv
 audit=resolved_source_groups.tsv annotations=accepted.annotations.tsv
 expectedinputs=10000 expectedreferences=21524 expectedmasked=2323
 baselineunmatched=1637 commoncutoff=RECORDED_MARGIN quantiles=64
 out=fresh_prefix t=1 -Xmx2g -ea -eoom

Fit three finite score-minus-own-GA thresholds with OR acceptance. Defaults are
refit=f upperidentity=.778: fit the 0.60 to below0.778 identity band and bypass
at0.778+. Missing and measured
impossible scores never pass a finite OR threshold; all missing rejects. The
separate common-cutoff comparator retains its own-CM/native-identity fallback.

The optional eventassignments=t export from analyzecobalamincmproposals.sh binds
native ordinal/worker/candidate IDs to the existing normalized one-to-one grade.
This tool requires all three target scores for every band event, including
unselected candidates, but fits only frozen emitted unmasked SELECTION calls.
By default, judging and masks never define quantiles or select the cutoff triple.

All-eligible saved-data refit (mode=fit only):
 refit=t upperidentity=.75
Pools both original SELECTION and JUDGING partitions for the quantile axes and
the same normalized objective. Masked calls still never fit; unavailable and
impossible scores still fail finite OR thresholds. Actual source-group IDs and
partition reports are retained, explicitly labelled not held-out evaluation.
All-eligible means all eligible emitted assignments in the supplied saved roster,
not a claim that unsupplied genomes or unselected raw events were fitted.
upperidentity accepts .75 or .778; omission keeps .778, even with refit=t.
The saved capture MUST still contain all three target scores for every event
below .778, including unselected events and the now-bypassed .75-to-.778 tail.
Those tail scores are validated but cannot define .75-fit quantiles. No new
alignment, CM scoring, capture or caller is performed by this mode.

The finite grid uses up to64 fitted-population quantiles per CM plus the recorded
common cutoff and a finite reject-all-fit-population control. Three-dimensional rejection
prefixes evaluate the complete tested surface. This is a grid approximation,
not an exhaustive continuous optimum. Ties use maximum minimum normalized
matched gain/unmatched reduction, then closest balance, then ascending cutoff
triple order. Objective scales are full-panel references/baseline unmatched.

Writes .axes.tsv, .surface.tsv.gz, .choice.tsv, .policy.tsv, .comparison.tsv,
.references.tsv, .genomes.tsv and .breakdown.tsv. The policy is header-only
unless both fitted-population measures improve. Every result describes fixed emitted
assignments, not a new caller path or a production recommendation. No rematching,
CM scoring, endpoint changes, reference relabeling or caller is performed.

Optional fixed comparison (mode=fit remains the default):
 mode=fixed fixedpolicy=accepted.result.policy.tsv
 gains=actual_original_or.comparison.loci.tsv expectedgains=283
 plus the same complete in/in2/raw/baseline/capture/sources/audit/annotations
 and population arguments above. commoncutoff and quantiles are not used.
Fixed mode rejects refit=t and upperidentity=.75. Its frozen .60-to-.778 score
coverage and target CM cutoffs are unchanged;
no fit, score pass or caller is run. Named rules are FULL_OR, PRUNED, A, B, C,
D1, D2, E and E_MINOR_66. A raises the cross floor to .63 main/.66 minor;
B lowers bypass to .75; C combines both. D1/D2 use native per-member
bypass/lower floors. E uses common .63 fallback and all three
CM tries; E_MINOR_66 uses native cross floors. E variants bypass at .75.
FULL_OR keeps all directions/native order; E variants keep all directions with
second/third/main priority; PRUNED/A/B/C/D1/D2 use direction mask30.
These fixed comparisons do not implement alignment-gated fallback F.

Fixed mode writes .fixed_policies.tsv, .fixed_comparison.tsv,
.fixed_references.tsv.gz, .fixed_genomes.tsv.gz, .fixed_events.tsv.gz and
.fixed_identity.tsv. Predictions include newly admitted capture-selected calls,
not only prior OR survivors. Unselected raw events remain separately counted;
no grade or reference assignment is invented for them. Original gained-reference
retention and captured member/model agreement are separate. Identity-tail tables
report matched/(matched+unmatched), masked counts and matched coverage separately
at .75/.778 for selection/judging/all. All grades use frozen assignments and do
not predict changed caller path selection.
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
 -cp "${DIR}current/" prok.SelectCobalaminClanAcceptance "$@"
