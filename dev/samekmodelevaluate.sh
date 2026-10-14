#!/bin/bash
# Exact candidate/control precision-recall curves from a frozen corpus fold.
# in=validation.tsv rows=master/rows.tsv net=model.bbnet fold=validation out=NEW_directory
# simdff=f default. Score >= cutoff accepts; targets are not inferred or resampled.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKModelEvaluate "$@"
