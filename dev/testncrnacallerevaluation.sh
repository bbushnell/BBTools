#!/usr/bin/env bash
# Finite native caller/lifecycle checks; optional child-java probe preserves the real JVM.
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
NCRNA_CALLER_REAL_JAVA=$(command -v java)
export NCRNA_CALLER_REAL_JAVA
if [[ -n "${NCRNA_CALLER_PROBE_BIN:-}" ]]; then
	PROBE_PATH="$NCRNA_CALLER_PROBE_BIN:$PATH"
	export PATH="$PROBE_PATH"
fi
exec "$NCRNA_CALLER_REAL_JAVA" $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.NcrnaCallerEvaluationTest "$@"
