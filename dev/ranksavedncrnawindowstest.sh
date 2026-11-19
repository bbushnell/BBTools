#!/bin/bash
# Independent tiny saved cache; no real genome or verifier work.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ranksavedncrnawindowstest.sh out=fresh_fixture_directory t=1 -Xmx512m -eoom
Checks literal two-model ranks, event-to-locus minima, panel-local identities,
zero-target input conservation, native-zero/no-window/no-compatible cases and
ambiguity against another scored subtype. Compares one versus two workers.
Rejects missing/duplicate sequences and methods, changed input/event/strand
bindings, wrong oracle model, false uniqueness and altered canonical tie tuples.
Synthetic correctness only; not an accuracy, timing or adoption result.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.RankSavedNcrnaWindowsTest "$@"
