#!/usr/bin/env bash
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CP=$DIR/current/
JAVA_BIN=$(command -v java)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
exec "$JAVA_BIN" $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" ml.CellNetSparseEncodingTest "$@"
