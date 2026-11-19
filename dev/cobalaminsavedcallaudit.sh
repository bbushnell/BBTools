#!/bin/bash
# Saved final-call correspondence, never a new caller or grader.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/cobalaminsavedcallaudit.sh in1=OFF_prefix in2=ON_prefix expectedinputs=1000 out=fresh_directory t=1 -Xmx1g -eoom
Reads complete normalized .genomes.tsv and .queries.tsv files, conserving every
declared matched/masked/unmatched call. Exact identity is input/contig/strand/
coordinates; member/model changes do not create a new exact call.
Unmatched ON calls partition into exact retained unmatched, no OFF overlap,
unique one-to-one OFF masked/matched overlap, or explicit other categories
(prior unmatched, multiple/mixed counterparts, exact-coordinate grade change).
Physical overlap is at least one base on the same input, contig and strand.
It is descriptive evidence, not a proved causal origin or a new truth match.
members.tsv reconciles additions/departures and label switches; categories.tsv,
on_calls.tsv, on_unmatched.tsv, off_calls.tsv and overlaps.tsv retain witnesses.
complete.tsv is written last. No sequence, model, caller or grader is executed.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$DIR/current/" prok.CobalaminSavedCallAudit "$@"
