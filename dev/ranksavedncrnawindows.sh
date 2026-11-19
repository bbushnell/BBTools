#!/bin/bash
# Saved-window ranking only; never reads genomes or runs a verifier.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ranksavedncrnawindows.sh panels=panels.tsv selection=report.inputs.tsv
  members=members.fa ledger=routing.tsv assignments=assignments.tsv
  census=membership.models.tsv summary=membership.summary.tsv ref=consensus.fa
  out=fresh_prefix expectedgroups=3153 expectedmodels=508 t=16 -Xmx8g -eoom

Compares forward k=7..11, plain/middle-masked, 256/512/1024/2048 buckets,
SSU and raw-key scores. Canonical plain k=9/11 controls must reproduce every
saved event's rank, tie span and model before results are written.
Reads the selected saved cache and frozen construction members, never genomes.
Keeps not-seen, seen and unknown strata separate. Outputs methods, summary,
per-locus minima, orientation comparisons, last reachable-locus losses at50
and work accounting. Saved ranks are not caller detections or adoption results.
The complete output prefix must be fresh; overwrite/append/sampling are refused.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=8g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.RankSavedNcrnaWindows "$@"
