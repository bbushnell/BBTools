#!/bin/bash
# Positive whole-assembly sequence identity for an explicit candidate-pair table.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_assembly_compare.sh in=candidate_pairs.tsv.gz out=fresh_prefix expectedpairs=360 maxassemblybases=100m t=1 bgzip=t pigz=f zl=9 -Xmx2g -eoom
Reads only the FASTAs named by the complete, bounded candidate-pair table from
ncrna_construction_assemblies.sh. Requires their audited base/contig dimensions.
Exact comparison ignores contig names, order, letter case and contig orientation.
Preserves contig multiplicity and DNA IUPAC ambiguity symbols; U/gaps reject.
Uses complete byte equality, not a hash-only or short-FMN identity inference.
Nonmatches remain biologically unresolved; circular rotation, alternate assembly
fragmentation and assembly-version relationships are outside this criterion.
At most two complete assemblies are retained; maxassemblybases bounds each one.
Writes .pairs.tsv.gz and .summary.tsv. Existing outputs, duplicate pairs, changed
input bindings, lossy base conversion, unexpected counts and dimensions fail.
This is positive-identity evidence, not a frozen training split or negative proof.
Compile on Dori; no annotation scan, model training or resource replacement.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaConstructionAssemblyCompare "$@"
