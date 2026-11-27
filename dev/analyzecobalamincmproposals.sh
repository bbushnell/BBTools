#!/usr/bin/env bash
# Saved per-call CM evidence and exact normalized-grade joins; no sequence work.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/analyzecobalamincmproposals.sh in=trace1[,trace2] in2=grade1[,grade2] \
  names=contact[,fixed32] profiles=cpu1.cm_endpoints.tsv[,cpu2.cm_endpoints.tsv] \
  bindings=panel.sources.tsv out=fresh_prefix expectedinputs=1000
Requires version-4 CM proposal traces and normalized grades of those exact GFFs.
Selected calls join by exact source/contig/strand/model/final bounds.
Contig terminals conserve all candidate records, including unselected proposals;
independent CM profile totals also detect omission of a whole unselected contig.
The endpoint oracle is a unique scored reference eligible for the ORIGINAL post-trim span;
unassociated/ambiguous cases remain explicit. No grade or caller is replayed.
Writes .calls.tsv, .bins.tsv, .tests.tsv and .summary.tsv. Guard discrimination
uses only finite proposals already allowed by existing constraints, with a fixed
original oracle. Exact-end costs accompany Youden separation; this conditional
proposal analysis does not predict downstream guarded detection or path changes.

members=RF00174 restricts trace/profile and selected-query linkage to that actual
CM scope; all normalized grade queries still conserve the complete genome totals.
Omission retains the original requirement for all three cobalamin members.

mode=scoreguard writes only .score_guard_calls.tsv and .score_guard_summary.tsv.
It reports actual SCORE_REJECTED requests and reconciles their number against
the independent rejected_score profile field. The original interval is retained;
the finite below-GA proposal is scored against the SAME predetermined reference.
All/selected/unselected request scopes use the unique original-span oracle;
a separate selected-recorded-TP scope uses the saved grader assignment. Missing,
ambiguous, masked and unselected cases remain visible. Within-three includes
exact ends. Recorded TP assignments may be strict or inclusive positive-overlap
matches; they do not inherit the separate original-oracle half-overlap cutoff.
The stored reference is retained without rematching. The inexact proposal
classes distinguish both-exact, within-three but
not both-exact, and at least one end beyond three. Beyond three is not itself
proof of worsening. Accurate rejected proposals are counterfactual endpoints,
not emitted gains or predictions of a new path. No caller or alignment runs.

mode=acceptance audit=resolved_groups.tsv rawgrades=grade.overlap.raw.tsv writes .acceptance_calls.tsv,
.acceptance_curves.tsv, .acceptance_coverage.tsv and .acceptance_native.tsv.
It describes emitted calls
using the accepted verifier identity, CM bits minus GA, and CM/GA. Curves use
every selection-observed inclusive cutoff, pooling ties, separately for the traced family
scope and each consensus model. Source-group parity zero selects; parity one
judges. This mode selects no cutoff and performs no three-zone policy fit.
Matched-inexact calls are the grader's recorded residual TP_INEXACT subset of
matches, joined from that exact raw grade without rematching. A separate
boundary_different column records differing final/reference endpoints; a strict
TP can have differing endpoints. Masked calls remain separate.
Missing scores remain NA and missing-count columns, never zero or a rejection.

acceptancescores=events.tsv acceptancecosts=costs.tsv optionally supply the
all-three-member capture from TimeNcrnaConsensusPolicy. Supply one of each per
trace. Exact input/worker/candidate IDs bind scores to the native effective
interval and verifier identity; main endpoint reuse must preserve score bits.
Per-input member costs conserve all candidate events, including unselected ones.
The endpoint profile includes extra main score-only requests; minor score-only
workers are accounted separately. Capture files must retain every admitted
candidate, so filtered-caller tables cannot be used as calibration captures.

nativecm=t (default for acceptance) reads canonical_new_gff from bindings and
joins the SAME recorded reference ID and measured CM identity. Protected model
assignments never substitute for cm_models. Native bits are printed tenth-bit
scores on the reference span; our anchored score uses the caller's own span
and band. Their difference is descriptive, not a scorer-parity test. Overlay-only
references, missing attributes/hits and absent same-CM measurements stay explicit.
The native table reports mean, sample standard deviation, nearest-rank quantiles
and GA mismatches, separately for each split and its matched-inexact subset.
nativecm=f skips that optional evidence. HBM is explicitly unavailable in this
trace format; this mode does not rescore sequences to fabricate an HBM curve.
Other members present only in normalized grades remain in coverage totals.
These fixed-output filters predict neither new caller paths nor recall below
the identity gate of the measured caller. No sequence search runs.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.AnalyzeCobalaminCmProposals "$@"
