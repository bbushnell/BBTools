#!/usr/bin/env bash
# Offline current-package loader controls; optional family=rnasep|tmrna|srp_small|srp_large.
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.NcrnaSharedBoundaryLoaderTest "$@"
