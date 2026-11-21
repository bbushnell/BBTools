#!/usr/bin/env bash

set -eo pipefail

DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "$DIR"/javasetup.sh
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment

java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.CountGCBaselineTest "$@"
