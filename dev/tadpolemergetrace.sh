#!/bin/bash
# Development-only Tadpole launcher: graph merge diagnostics, no per-kmer logs.
# Accepts ordinary tadpole.sh arguments. All assembly decisions are unchanged.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="$DIR/current/"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" assemble.TadpoleMergeTrace "$@"
