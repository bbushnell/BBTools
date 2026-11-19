#!/bin/bash
# Native saved-window acceptance development tool.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/replay_ncrna_windows_test.sh out=fresh_fixture_directory -Xmx512m -eoom
Runs the real capture once, then checks new-context and legacy-context replay,
retained zero-window input, both strands, physical edges and failed control.
Correctness fixture only; fresh child JVMs isolate caller static configuration.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ];do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError \
	-cp "${DIR}current/" test.prok.ReplayNcrnaWindowsTest "$@"
