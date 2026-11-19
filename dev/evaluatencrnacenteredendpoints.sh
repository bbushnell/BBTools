#!/usr/bin/env bash
# Compare saved predictions in the original coordinate frame, retaining holds.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/evaluatencrnacenteredendpoints.sh in=baseline.scores.tsv \
 centered=centered.scores.tsv ledger=r10.evaluation.5prime.centered.tsv \
 loci=native.loci.tsv end=5prime expected=13401 usable=13000 out=fresh_prefix

Reuses saved 69-input/22-output scores from EvaluateNcrnaPositionalNet fullerrors=t.
No inference or fitting. Every original trial remains in the denominator.
usable= is the accepted number of VECTORIZED ledger rows. If zero, omit centered=;
all trials still appear with their explicit hold and baseline predictions.
Both score streams require their original provenance and terminal COMPLETE marker.
The native locus table supplies observed context lengths, NOT whole-contig lengths.

Outputs: .summary.tsv and .trials.tsv. Reports raw, old-network margins 0/0.1,
center-only and centered-network margins 0/0.1 in original coordinates. Invalid
centers revert to the raw pair; valid centers without features retain the center.
Final network proposals are scored without a new physical hold, matching the old
independent-end assay; crossed, outside-context and too-short proposals are counted
per arm. Caller joint holds remain a separate integration test. Counts are trials,
not unique genes or caller detections. Unsupported offsets and already-invalid
raw spans remain explicit. All outputs must be fresh.

Geometry-only diagnostic (no score inputs, no inference or metric replay):
  geometryonly=t ledger=all_trials.tsv loci=native.loci.tsv end=5prime \
  expected=13401 out=fresh_prefix
Writes .geometry.tsv with the raw-by-centered geometry cross-tab, distinguishing
already-invalid holds, newly invalid centers and repaired raw spans. Omit in=,
centered= and usable= in this mode. Context limits are not whole-contig lengths.

Finite controls: dev/evaluatencrnacenteredendpoints.sh test out=fresh_directory
HELP
 exit 0
fi
target=prok.EvaluateNcrnaCenteredEndpoints
if [[ $1 == test ]];then
 target=test.prok.EvaluateNcrnaCenteredEndpointsTest
 shift
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$target" "$@"
