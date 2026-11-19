#!/bin/bash
# Reconstruct only saved model orders; no genome reader or verification.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/hybridsavedncrnawindows.sh panels=panels.tsv selection=report.inputs.tsv
  grid=saved_grid_prefix method=64 members=members.fa ledger=routing.tsv
  assignments=assignments.tsv census=membership.models.tsv
  summary=membership.summary.tsv ref=consensus.fa out=fresh_prefix
  expectedgroups=3153 expectedmodels=508 t=16 -Xmx8g -eoom

Reads original native shared-key counts and reconstructs the selected sketch.
Proposes native[0], sketch[0], native[1], sketch[1], skipping repeated models;
duplicates consume proposal turns but never another attempted-model rank.
Native zero-support models are excluded from that side; the sketch is complete.
Checks exact native/canonical event rank/tie/model tuples and saved minima.
The selected new sketch has saved locus-minimum parity only: an independent
new-method event-tuple cache does not exist. This limitation is in the work table.
Outputs separate unseen/seen/unknown ranks, per-input reachability/losses and
method metadata. Saved ranks are not caller detections, speed or adoption.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.HybridSavedNcrnaWindows "$@"
