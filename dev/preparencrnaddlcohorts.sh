#!/bin/bash
# Freeze saved-window identities once, reusing accepted values when available.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/preparencrnaddlcohorts.sh panels=effective.panels.tsv selection=report.inputs.tsv ref=consensus.fa out=fresh_prefix expectedmodels=508 expectedgroups=0 t=1 -Xmx4g -eoom
Reads the original saved windows and native oracle. No caller or sketch fitting.
expectedgroups=0 counts the complete selected inventory before validating it.
Optional reuse requires previous=accepted_prefix previouspanels=old.panels.tsv
previousselection=old.inputs.tsv previousgroups=N together. The same frozen
consensus must be used. Prior values and missing dispositions are restored only
for identical source/input/parent/coordinate and native model/event witnesses.
Only new parents reach identity measurement. Writes .identity.tsv, .groups.tsv,
.identity_bins.tsv and .preparation.tsv with reused/new/missing counts.
The historical not-seen cohort labels are preserved; the replicated scorer
derives its full-comparable low-identity cutoff before scoring from these values.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ];do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
CP="${DIR}current/"
if [ -f "${DIR}bbtools.jar" ];then CP="${DIR}bbtools.jar";fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=4g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.PrepareNcrnaDdlCohorts "$@"
