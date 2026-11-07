#!/bin/bash
# Finite saved-hit projection tests. Pass work=a_fresh_owned_directory.
set -eo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
exec java -ea -XX:+ExitOnOutOfMemoryError $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.NcrnaCmHitProjectorTest "$@"
