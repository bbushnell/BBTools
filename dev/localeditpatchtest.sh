#!/bin/bash
# Run deterministic patch transaction, explicit-region and automatic-region tests.
# JVM options are parsed here; all three test mains receive no application arguments.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
for test in LocalEditPatchBatchTest LocalEditPatchCorrectorTest LocalEditPatchRegionsTest;do
  java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" "test.assemble.$test"
done
