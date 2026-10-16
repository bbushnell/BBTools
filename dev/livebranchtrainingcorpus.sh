#!/bin/bash
# Streams frozen12-input exports into TRAIN and separate VALID populations.
set -eo pipefail
if [[ $# -eq 0 ]]; then
    echo 'Usage: livebranchtrainingcorpus.sh in=selected.tsv expected=848 out=NEW_directory'
    echo 'TRAIN must use controlfraction=0.1 and its frozen genome seed; VALID must be unsampled.'
    echo 'Writes train.tsv, validation.tsv (grey only), validation_controls.tsv, rows.tsv, census.tsv and complete.txt.'
    exit 0
fi
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.LiveBranchTrainingCorpus "$@"
