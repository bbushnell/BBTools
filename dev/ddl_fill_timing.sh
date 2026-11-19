#!/bin/bash
# One-thread native versus shared-iterator fill timing; every sketch has 512 bins.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
    cat <<'HELP'
Usage: dev/ddl_fill_timing.sh panels=saved_panels.tsv selection=saved_inputs.tsv out=fresh.tsv seconds=10 t=1 -Xmx4g -eoom
Compares native canonical k7/k9 filling with forEachKmer + hashAndStore on
identical saved windows, after complete key parity. Three alternating rounds,
at least ten main-thread CPU seconds per arm. The actual forward/masked/tagged
dual fill is a separate workload. Includes clearing and consuming all keys;
excludes loading and construction. seconds=0 with no inputs runs finite parity.
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
parseJavaArgs --xmx=4g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.cardinality.DdlFillTiming "$@"
