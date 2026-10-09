#!/usr/bin/env bash
# Bounded synthetic two-pass refinement and graph compatibility fixtures.
set -eo pipefail
DIR=$(cd "$(dirname "$0")/.." && pwd)
CP=$DIR/current
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
 case "$arg" in -Xmx*|-Xms*|-ea|-da|-eoom) JVM_ARGS+=("$arg");; *) ARGS+=("$arg");; esac
done
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" prot.AAGraphTest
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" prot.HbmProfileRefinerTest "${ARGS[@]}"
