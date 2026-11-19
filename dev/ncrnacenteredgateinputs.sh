#!/usr/bin/env bash
# Combine real held-call35-input gate matrices; no fit or scaling.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/ncrnacenteredgateinputs.sh in=accepted_real_fold_manifest.tsv \
 member=RF00174 end=5prime out=fresh_directory -Xmx1g -eoom
Manifest columns: fold, dataset, prefix, expected_calls, expected_proposals (TAB).
Each prefix binds one seed58 real gate adapter's .rows.tsv.gz/.vectors.tsv.gz.
Folds0..3 become train.tsv.gz; fold4 becomes monitor.tsv.gz. Every fold must
contribute real proposals; outer rows, duplicate physical references, truncated
streams, shape/label/displacement disagreement fail. No transform is fitted;
RegressionTrainer later learns normalization from train.tsv.gz alone.
Outputs are BGZF matrices plus per-dataset and per-fold population accounting.
The monitor is explicitly not nested independent. No gate fit, caller or cutoff.
Finite controls: put test first and out=fresh_directory.
HELP
 exit 0
fi
target=prok.NcrnaCenteredGateInputs
if [[ $1 == test ]];then
 target=test.prok.NcrnaCenteredGateInputsTest
 shift
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$target" "$@"
