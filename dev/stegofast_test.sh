#!/bin/bash
set -e
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
. "$DIR/memdetect.sh"
parseJavaArgs "--xmx=512m" "--xms=128m" "--mode=fixed" "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.stego.StegoFastTest "$DIR"
