#!/bin/bash
# Literal saved-input and exact-sequence panel-binding fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_caller_panel_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Checks original and combined populations, independent source slots, zero-target
inputs, unrankable and UNKNOWN parents, UTF-8 paths and distinct sequence/source
states. Rejects missing/duplicate parents, wrong source slots and duplicate
assembly groups. All genome paths are nonexistent; metadata only, no caller.
Compile complete affected packages on Dori before running this fixture.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaDdlCallerPanelTest "$@"
