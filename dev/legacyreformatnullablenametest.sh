#!/usr/bin/env bash
# Direct legacy Read-overload nullable-ID controls; no driver constructors or CLI workers.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.LegacyReformatNullableNameTest "$@"
