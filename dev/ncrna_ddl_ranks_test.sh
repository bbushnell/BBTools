#!/bin/bash
# Literal deterministic ranks, complete score ties and native whole-window scores.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_ranks_test.sh t=1 -Xmx512m -eoom
Checks literal ranks/ties, absent oracle successes, native SSU float scores,
one-bucket matches, query replacement and borrowed input ownership.
Also checks forward/plain and forward/middle-masked words for k7–11 and widths
256/512/1024/2048: a sole central SNP scores zero plain and positive masked;
even k clears the right central base while preserving a left central SNP.
Canonical words exactly match the native sketch, including N/lower-case bases.
Synthetic correctness only; no genome accuracy or performance result.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaDdlRanksTest "$@"
