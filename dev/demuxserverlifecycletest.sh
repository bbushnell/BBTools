#!/bin/bash
# One disposable-local DemuxServer lifecycle case: mode=baseline|fixed case=<name>.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
if [[ "${1:-}" == "--child" ]]; then
	shift
	exec java -ea -XX:+ExitOnOutOfMemoryError -Xmx128m -Xms32m \
		-cp "$DIR/current/" test.barcode.DemuxServerLifecycleTest child=t "$@"
fi
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.DemuxServerLifecycleTest \
	installation="$DIR" "$@"
