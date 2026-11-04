#!/usr/bin/env bash
# Five finite real-public-tool controls; no resource corpus or genome evaluation.
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" \
	test.prok.NcrnaSequentialRegenerationTest "$@"
