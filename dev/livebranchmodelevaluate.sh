#!/bin/bash
# Exact livebranch12-v1 curves, with grey candidates and controls kept separate.
# in=validation.tsv rows=master/rows.tsv net=model.bbnet partition=validation out=NEW_directory
# partition=train|validation|validation_controls; simdff=f default; simdfma follows the standard parser.
# Use simdff=f simdfma=f for scalar inference. Accepts score>=cutoff; never resamples.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.LiveBranchModelEvaluate "$@"
