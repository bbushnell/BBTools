#!/bin/bash
# Launch the scalar BubblePopper retained-span probe.
# The probe verifies span arithmetic without constructing contigs.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=64m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.BubblePopperRetainedSpanProbe "$@"
