#!/usr/bin/env bash
# Finite saved-call/paired-library controls; no genome calling or resource build.
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" \
	test.prok.NcrnaCallerReversionTest "$@"
