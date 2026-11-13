#!/bin/bash
# Finite barcode filtering/quality checks; one case/scenario per child JVM.
# mode=baseline|fixed case=filter scenario=none|minimum|average|both|single dir=NEW
# mode=baseline|fixed case=quality scenario=control|barcode|read1|read2|high dir=NEW
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.CorrelateBarcodesRegressionTest "$@"
