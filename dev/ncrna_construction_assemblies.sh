#!/bin/bash
# Known assembly aliases and dimension-compatible comparison candidates.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_assemblies.sh genomes=original.genomes.tsv evalcensus=eval.genomes.tsv evalmembers=eval.memberships.tsv newcensus=new.genomes.tsv newmembers=new.memberships.tsv expectedgenomes=39979 expectedeval=11769 expectednew=29313 out=fresh_prefix t=1 bgzip=t pigz=f zl=9 -Xmx2g -eoom
Reads metadata and resolves existing FASTA paths without reading sequence payloads.
Original zero-locus rows have unmeasured zero dimensions; the positive census supplies
their genome sizes. Measured original dimensions must still match the census.
Canonical file aliases and explicit assembly-record identities define known groups.
Membership claims must bind to the census file by canonical identity; foreign or
unavailable paths remain unresolved. Retained paths are emitted as UTF-8 bytes.
Equal base/contig counts nominate whole-assembly comparisons, never prove identity.
Short FMN-sequence identity is not consulted and cannot exclude another genome.
Writes .inputs.tsv.gz, .candidate_pairs.tsv.gz and .summary.tsv. Existing outputs fail.
This is an assembly audit, not a frozen member split or sequence-equivalence proof.
Compile complete prok on Dori first; source files must exist at their recorded paths.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaConstructionAssemblies "$@"
