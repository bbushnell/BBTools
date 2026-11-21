#!/usr/bin/env bash
# SAM breaklength baseline/strand-repair controls; fresh JVM and outdir.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.ReformatSamBreakBaselineTest "$@"
