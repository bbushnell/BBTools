#!/bin/bash
# One bounded TagAndMerge failure case: mode=baseline|fixed case=<name>.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
if [[ "${1:-}" == "--tag-child" ]]; then
	shift
	exec java -ea -XX:+ExitOnOutOfMemoryError -Xmx128m -Xms32m \
		-cp "$DIR/current/" barcode.TagAndMerge "$@"
fi
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.TagAndMergeFailureTest \
	installation="$DIR" "$@"
