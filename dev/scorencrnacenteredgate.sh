#!/usr/bin/env bash
# Cache one inference per real whole-move proposal; no cutoff selection.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/scorencrnacenteredgate.sh in=real_gate_manifest.tsv net=gate.bbnet \
 member=RF00174 end=5prime seed=58 partition=internal out=fresh_directory
The manifest has fold,dataset,prefix,expected_calls,expected_proposals columns.
Each prefix names a complete real adapter matrix, all-call audit and first-stage
score cache. Internal folds use auxiliary first-stage seed58; outer folds=-1
require the matching full first-stage seed (58/59/60) for the gate seed.
Output scores.tsv.gz keeps every original call, all raw first-stage score bits
and one gate score per proposal. Holds never invoke the gate. No fit, first-stage
inference, cutoff selection, caller or matching is performed. BGZF output.
Finite controls: put test first and out=fresh_directory.
HELP
 exit 0
fi
target=prok.ScoreNcrnaCenteredGate
if [[ $1 == test ]];then
 target=test.prok.ScoreNcrnaCenteredGateTest
 shift
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$target" "$@"
