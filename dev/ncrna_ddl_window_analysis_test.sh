#!/bin/bash
# Synchronous ordinary-window analysis, literal spans and observer parity.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_window_analysis_test.sh t=1 -Xmx512m -eoom
Checks plus/minus native seed inputs, full-oracle literal spans, both profile
geometries, complete ordinary GFF/attempt parity, unchanged borrowed bytes,
empty contigs, repeated event coordinates and duplicate native seqid rejection.
The native ncRNA CPU profiler must be disabled; analysis owns separate counters.
Synthetic correctness only; no biological ranking or speedup claim.
Compile complete prok/test.prok packages on Dori first.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaDdlWindowAnalysisTest "$@"
