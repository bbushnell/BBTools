#!/usr/bin/env bash
# Conditional magnitude-only replay of source-qualified saved real calls.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/ncrnaoffsetthresholdreplay.sh in=geometry_datasets.tsv \
 folds=complete_assemblies.tsv ref=frozen_consensus.fa member=RF00174 \
 out=fresh_directory -Xmx2g -eoom

Dataset TSV header: dataset, role, geometry, expected_rows (tabs).
Roles are TRAIN or OUTER; supply the complete support pool and two outer panels.
Geometry caches come from evaluatencrnarealendpoints.sh mode=geometry.
Complete assembly folds come from NcrnaEndpointFoldInputs allassemblies=t.
Every support/evaluation assembly must map explicitly. All outer assemblies
are excluded from learning. Each internal fold's offsets exclude that fold.

Support floor10, nearest-rank lower medians. Apply a supported median only when
abs(median)>T. One threshold controls both ends; fixed ends hold. Invalid joint
spans restore raw bounds. Single-end scopes hold the other endpoint raw.
Sweep distinct fold-trained median magnitudes,0,and never-shift. Fold4 alone
selects by greatest end better-minus-worse, then fewer worse ends, then fewer
moves; retain all ties. Apply frozen choices to outer panels with full-training
offsets. Outer truth never selects thresholds. No networks or caller run.

Outputs: five held-fold and full offset resources, threshold choices, complete
internal curves, outer chosen/control rows and population accounting. Report
full and original-near RMS, exact/within3/10, both ends, B/W/T and physical holds.
This is conditional matched-call ability, not detection, path reproduction or
nested independent monitoring. Finite controls: put test first and out=fresh_dir.

Gate preprocessing: mode=residuals offsets=accepted_threshold_result_directory
reuses the five held-fold and full offset resources without fitting medians or
repeating selection. It reads TRAIN geometry only (never OUTER files), checks
model order/support against the same fold population, and exports signed
residual histograms with actual physical/fixed/minimum-length joint guards.
Use the same in=, folds=, ref= and member= bindings and a fresh out= directory.
HELP
 exit 0
fi
target=prok.NcrnaOffsetThresholdReplay
if [[ $1 == test ]];then
 target=test.prok.NcrnaOffsetThresholdReplayTest
 shift
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$target" "$@"
