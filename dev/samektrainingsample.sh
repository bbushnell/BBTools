#!/bin/bash
# Preserves all known training candidates and uniformly samples an exact control quota.
# in=master/train.tsv rows=master/rows.tsv out=NEW_directory controls=N seed=N
# Writes train.tsv, rows.tsv (with original row numbers), complete.txt.
# Never opens validation/test tables. Input must be a validated SameKTrainingCorpus.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKTrainingSample "$@"
