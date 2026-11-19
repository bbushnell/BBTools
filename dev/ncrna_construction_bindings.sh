#!/bin/bash
# Original FMN construction identities; no eligibility, extraction or training.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_bindings.sh genomes=old.genomes.tsv loci=old.loci.tsv.gz assignments=old.assignments.tsv evalcensus=eval.genomes.tsv evalmembers=eval.memberships.tsv newcensus=new.genomes.tsv newmembers=new.memberships.tsv expectedgenomes=39979 expectedmembers=26896 out=fresh_prefix t=1 -Xmx2g -eoom
Reads existing metadata only. Preserves all old seed inputs, their original model
assignments and every literal cohort/input binding. Assembly claims are explicit;
sequence-alias status remains UNVERIFIED and eligibility remains NOT_EVALUATED.
Writes .members.tsv.gz, .bindings.tsv and .summary.tsv. Existing outputs fail.
Compile the complete prok package on Dori first. No source payload is regenerated.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaConstructionBindings "$@"
