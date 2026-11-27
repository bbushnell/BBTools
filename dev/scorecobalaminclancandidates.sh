#!/usr/bin/env bash
# Saved-candidate covariance scoring; developer measurement, no caller run.
set -eo pipefail
if [[ $# == 0 || ${1:-} == -h || ${1:-} == --help ]];then
 cat <<'HELP'
Score every saved cobalamin candidate in the explicit identity band with all
three covariance models. Reuse its own-CM measurement; select each other CM's
anchor from that member's existing 1024-bin profiles (best bucket count, first
library ordinal on ties), then reconstruct an independent Quantum trace.
The DDL query is the saved post-refinement candidate span, sense-oriented.
Main uses its existing learned corridor/input padding; minors use radius32,
retry64 on band contact, and20-base flanks. No endpoints or path are changed.

Usage: dev/scorecobalaminclancandidates.sh in=work.acceptance.tsv \
 in2=work.inputs.tsv familyroot=selected_group scoreroot=score_group \
 profiles=main.tsv,second.tsv,third.tsv expectedinputs=10000 \
 lowidentity=.60 highidentity=.778 out=fresh_prefix t=16 -Xmx4g

The launcher enables the native CPU clock; the output clock property is only
a diagnostic switch here, not an additional output file. Writes .scores.tsv
(three rows per band candidate, retaining its original event and index) and
.costs.tsv (emitted-member by scoring-member counts and added CPU nanoseconds).
Own-score reuse adds zero CPU to the recorded comparison. Missing maps/anchors
remain explicit unavailable scores; impossible scores remain negative infinity.
These are candidate events, not selected loci or a new reference grade. Row
order across workers is unspecified; event index plus scoring member is stable.
No scoring cutoff is fitted, installed or interpreted by this program.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -Dbbtools.diagnostic.ncrnaCpuOut=ENABLED_ONLY_NO_CALLER \
 -cp "$DIR/current/" prok.ScoreCobalaminClanCandidates "$@"
