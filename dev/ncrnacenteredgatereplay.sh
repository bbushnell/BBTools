#!/usr/bin/env bash
# Physical single/joint coordinate replay from frozen real gate scores.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/ncrnacenteredgatereplay.sh mode=select in=start/scores.tsv.gz \
 in2=stop/scores.tsv.gz member=RF00174 seed=58 expected=full_internal_call_count \
 out=fresh_directory -Xmx1g -eoom
mode=evaluate additionally requires choices=accepted_monitor/choices.tsv and
two OUTER caches from one original panel with its exact expected=call_count.

Selection uses fold4 only. Per-end score>=cutoff maximizes better-minus-worse,
then fewer harms, then fewer actual moves with the other end RAW and physical
guards. Every tied optimum remains. Matched-rate controls retain tied-score
blocks and explicit unavailable brackets. Never retune outer movement counts.

summary.tsv keeps raw, center-only, centered margins0/.1, all whole proposals,
diagnostic ignored-head and frozen gate policies. Single-other-RAW and paired
joint outcomes remain separate. Invalid gate pairs revert RAW; existing neural
control pairs revert their guarded center. All matched calls remain counted.
Full/raw-near RMS, signed changes, both-end counts, harms and physical holds are
reported with raw sums for later member/clan accounting. No inference or fit.
This is conditional endpoint ability, not caller/path reproduction or adoption.
Finite controls: put test first with out=fresh_directory.
HELP
 exit 0
fi
target=prok.NcrnaCenteredGateReplay
if [[ $1 == test ]];then
 target=test.prok.NcrnaCenteredGateReplayTest
 shift
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$target" "$@"
