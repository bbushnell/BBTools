#!/bin/bash
# Run deterministic graph-only unit tests; optional argument writes CLI fixtures.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="$DIR/current/"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" assemble.TadpoleGraphTest "$@"
for test in BubblePopperUnitTest PathPreservingBubbleSimplifierSpec ReadThreadedXResolverUnitTest \
    CrossKTipOverlapperUnitTest SimpleOmnitigExtractorUnitTest TadpoleMultiUnitTest ContigGraphClassifierUnitTest; do
  java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" "assemble.$test"
done
