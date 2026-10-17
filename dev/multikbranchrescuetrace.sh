#!/bin/bash
# Development-only multi-K rescue merge tracing; accepts multikbranchrescue.sh flags.
# Existing graph/overlap diagnostics only, without per-kmer traversal logging.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=4g --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.MultiKBranchRescueTrace "$@"
