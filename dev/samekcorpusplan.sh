#!/bin/bash
# Freezes group-disjoint 70/15/15 splits (rounded by group) and one random K
# per genome from a Plan35 samples.tsv. No training or sequence inspection.
# Usage: in=samples.tsv out=new.tsv root=/absolute/data expected=1000
# splitseed=370194 kseed=370195. Outputs must be new. Existing input files
# are NOT verified by this planner; validate them before running collection.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKCorpusPlan "$@"
