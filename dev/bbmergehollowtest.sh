#!/bin/bash
# Exercises BBMerge constructor routing and its real paired extension helper.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.BBMergeHollowTest in="$DIR/testdata/bbmerge_hollow_pairs.fa" "$@"
