#!/bin/bash
# Literal saved-rank sums, named cases and explicit seen/unseen/unknown strata.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_saved_rank_summary_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Checks literal common-locus sums, means, coverage/quantiles, seen/unseen/unknown
denominators, named zero-support/unresolved cases, UTF-8 paths and CLI parity.
Truncated, misordered, misbound, missing-rank and wrong-count reports fail loudly.
Synthetic metadata only; no sequences, caller, alignment or hybrid reconstruction.
Compile the full canonical current/ tree on Dori before running.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaSavedRankSummaryTest "$@"
