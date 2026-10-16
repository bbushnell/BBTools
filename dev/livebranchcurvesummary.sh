#!/bin/bash
# Selects maximal true rescues under fixed false-rescue counts from a complete exact curve.
# in=grey.tsv out=NEW.tsv expected=42108; ties keep the higher cutoff. No TEST access.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.LiveBranchCurveSummary "$@"
