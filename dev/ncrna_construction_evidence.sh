#!/bin/bash
# Saved FMN construction evidence; no eligibility, extraction or training.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_evidence.sh bindings=old_prefix evaluation=eval_tiers discovery=new_tiers expectedmembers=26896 out=fresh_prefix t=1 bgzip=t pigz=f zl=9 -Xmx2g -eoom
Reads existing bindings and native source/event correspondence. Preserves original
member IDs, cohort namespaces, all FMN events and measured disjoint-core overlap.
Source-line and direct bound-coordinate links are distinguished explicitly;
missing saved GFF witnesses remain flagged even when coordinate evidence exists.
Other-family support is diagnostic; eligibility remains NOT_EVALUATED.
Writes .old_evidence.tsv.gz, .links.tsv.gz, .fmn_events.tsv.gz and .summary.tsv.
Existing outputs fail. No sequence-alias proof or biological admission is performed.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaConstructionEvidence "$@"
