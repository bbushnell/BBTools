#!/bin/bash
# Subprocess fixture: expected exit1 with original worker OOM and memory guidance.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	echo 'Usage: dev/callgenes_worker_oom_test.sh in=tiny.fa t=1 -Xmx32m'
	echo 'Do not use -eoom; the ordinary worker catch,not JVM startup,must handle OOM.'
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
CP="${DIR}current/"
if [ -f "${DIR}bbtools.jar" ]; then CP="${DIR}bbtools.jar"; fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=32m --xms=16m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" test.prok.CallGenesWorkerOomTest "$@"
