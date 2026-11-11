#!/bin/bash
# Run deterministic graph-only unit tests; optional argument writes CLI fixtures.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="$DIR/current/"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" test.assemble.TadpoleGraphTest "$@"
for test in test.assemble.BubblePopperUnitTest test.assemble.PathPreservingBubbleSimplifierSpec \
    test.assemble.ReadThreadedXResolverUnitTest test.assemble.CrossKTipOverlapperUnitTest; do
  java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" "$test"
done
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" test.assemble.FusionKmerSupportTest
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" test.assemble.SimpleOmnitigExtractorUnitTest
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" test.assemble.TadpoleMultiUnitTest
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" test.assemble.ContigGraphClassifierUnitTest
