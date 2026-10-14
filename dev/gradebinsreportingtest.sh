#!/usr/bin/env bash
# Shared GradeBins summary regression, including QuickBin-facing methods.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
  case ${arg,,} in
    -ea|-da|-eoom|eoom|-xmx*|-xms*|simd|simd=t|simd=f) JVM_ARGS+=("$arg");;
    *) ARGS+=("$arg");;
  esac
done
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$DIR/current" bin.GradeBinsReportingTest "${ARGS[@]}"
