#!/usr/bin/env bash
# Constructor-only RenameReads observation; fresh JVM and disposable inputs.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.RenameReadsConstructorTest "$@"
