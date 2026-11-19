#!/bin/bash
# Integrate saved positive assembly matches without reopening genome payloads.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_assembly_groups.sh in=audit.inputs.tsv.gz nominations=audit.candidate_pairs.tsv.gz comparisons=comparison.pairs.tsv.gz out=fresh_prefix expectedeval=11769 expectednew=29313 expectedpairs=360 t=1 bgzip=t pigz=f zl=9 -Xmx2g -eoom
Reads saved metadata only. Preserves every original input and comparison outcome.
Each comparison must retain its exact nomination row, order, IDs, paths and criterion.
Only positive exact matches join existing assembly groups. Nonmatches never join
groups or prove biological independence; an independent positive match may still
exclude that input. Unresolved pair flags remain visible even in excluded groups.
Reports unique discovery dispositions separately from pair counts. No member split
is frozen and no short-gene identity, extraction, training or genome reread occurs.
Writes .inputs.tsv.gz, .pairs.tsv.gz and .summary.tsv; existing outputs fail.
Compile full current/ on Dori first.
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
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaConstructionAssemblyGroups "$@"
