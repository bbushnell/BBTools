#!/bin/bash
# Streams samek27-v1 collection files into native #dims 27 1 trainer tables.
# in=cohort.tsv inventory=inventory.tsv root=collection_root out=NEW_directory
# expected=1000. Unknown labels are excluded and counted; folds are never redrawn.
# Output: train/validation/test.tsv, rows.tsv, census.tsv, complete.txt.
# Only a completed directory is usable: an input failure may leave partial files
# but never complete.txt. No model is fitted and no source files are changed.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKTrainingCorpus "$@"
