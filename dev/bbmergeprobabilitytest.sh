#!/usr/bin/env bash
# Finite probability checks: phase=baseline outdir=NEW or
# phase=fixed outdir=NEW reference=BASELINE_FIXTURE_DIRECTORY (the runner's cases/).
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea simd=f "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.BBMergeProbabilityTest "$@"
