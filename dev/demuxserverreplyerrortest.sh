#!/bin/bash
# In-memory DemuxServer reply writer failure fixture: mode=baseline|fixed case=<name>.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.DemuxServerReplyErrorTest "$@"
