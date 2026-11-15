#!/usr/bin/env bash
# Finite hiseq equals-in-file-argument parser checks; constructs tools but never processes reads.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed -ea -eoom simd=f
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$DIR/current" \
 test.hiseq.HiseqEqualsPathFixtureTest "$@"
