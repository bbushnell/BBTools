#!/usr/bin/env bash
# Offline comparison harness; no production assignment or alignment algorithm changes.
set -eo pipefail
kind=${1:?data or run}; shift
case $kind in data) main=prot.HmmComparisonData;; run) main=prot.HmmComparisonRunner;; *) exit 2;; esac
DIR=$(cd "$(dirname "$0")/.." && pwd)
CP=$DIR/current
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
	case "$arg" in -Xmx*|-Xms*|-ea|-da|-eoom) JVM_ARGS+=("$arg");; *) ARGS+=("$arg");; esac
done
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$main" "${ARGS[@]}"
