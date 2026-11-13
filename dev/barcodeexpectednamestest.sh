#!/bin/bash
# Finite expected-name cases: mode=baseline|fixed case=annotated|plain|rc|dual|tab|generic source=file|literal dir=NEW.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.BarcodeExpectedNamesTest "$@"
