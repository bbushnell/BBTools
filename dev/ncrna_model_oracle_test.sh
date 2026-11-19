#!/bin/bash
# Finite single-model verification through the ordinary RF00174 verifier.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_model_oracle_test.sh -Xmx512m -eoom
Checks original model identity, interior and physically clipped coordinates,
both caller-oriented strands, recomputed seed positions, HBM rescue and
independence from a rejecting mapping index. Requires assertions enabled.
Also checks cached exhaustive results, original-forward spans, described headers,
exact input identity, replacement and ambiguous tandem compatibility links.
No genome input, learned-resource result or full ranking experiment.
Compile complete prok and test.prok packages on Dori before running.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaModelOracleTest "$@"
