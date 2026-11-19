#!/usr/bin/env bash
# Diagnostic only: raw state/output bits; no supplied network is trained/exported.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Capture exact network state and predictions for pre/post-merge comparison.
Usage: capturencrnamlcompatibility.sh in=networks.tsv out=fresh_capture.tsv rows=32 ow=f
Manifest: id,net,frame,vectors,inputs,outputs,provenance,seed (tabs).
frame=-1 selects a standalone file; 0 and above select explicit bundle frames.
REAL/-1 rows use saved raw matrices of exactly the declared input width.
SYNTHETIC/nonnegative-seed rows use a named source-bound fixture:
synthetic:trna_profiles_v1, synthetic:rrna_planted_r12_v1 or
synthetic:rrna_planted_r16_v1. Generated rows and extremes are explicitly labeled;
their sampling intervals are fixture choices, not universal biological bounds.
Records ordinary, worker-copy, inference-only and dense-optimized modes separately,
plus deterministic input edges and one small in-memory seeded tanh update.
legacy=t keeps that default training control and requires saved ten-input rows.
legacy=f omits only that control for a network-only diagnostic of any input width.
No scientific training, network export, selection or adoption is performed.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.NcrnaMlCompatibilityCapture "$@"
