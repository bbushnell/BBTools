#!/bin/bash
# Ordinary single-member seed/root loading and shared-index regression fixture.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	echo 'Usage: dev/callgenes_cobalamin_seed_override_test.sh out=fresh_fixture_directory member=cobalamin_rf00174 root=t|f t=1 -Xmx2g -eoom'
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
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" test.prok.CallGenesCobalaminSeedOverrideTest "$@"
