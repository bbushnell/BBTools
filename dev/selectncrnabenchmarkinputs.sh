#!/bin/bash
# Family-blind sampling of saved whole-input census metadata.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden

Usage: dev/selectncrnabenchmarkinputs.sh in=genomes.tsv memberships=memberships.tsv \
       taxonomy=taxonomy.tsv sizes=1000,11769 out=new_prefix seed=20261111

Selects benchmark inputs from a saved census; reads no genome, annotation or
reference payloads. Selection uses only domain/phylum, nucleotide bases and
FASTA-record counts. The taxonomy must bind every exact census input/path once.
Unknown taxonomy remains an explicit stratum. No family-positive counts are used.

Each domain/phylum receives one slot. Remaining slots are proportional to its
remaining population, with deterministic largest-remainder allocation. Within
each group, systematic sampling over joint log2(base-count)/log2(record-count)
bins preserves workload variety; seeded ranks order members within a bin.
The minimum-per-group rule intentionally overrepresents small groups. Inspect
the composition table before describing the result as representative.

sizes=       Exact positive sizes, each between the number of groups and census
             population. A full-population size keeps every original input.
seed=20261111  Retained ranks use Java's seeded Random in original census order.

For each size N, writes prefix.N.inputs.txt, .inputs.config (one comma-list),
.flat.config (one in= record per input), .selection.tsv, .composition.tsv,
and .totals.tsv (corpus/selected input, nucleotide-base and FASTA-record totals).
Output order is the original census order. All outputs must be fresh; ow=f.
Input paths may contain UTF-8 and equals signs, but not comma/record delimiters.
Zero exit and the final completion marker are required; failures may leave
partial outputs. Hash bindings and corpus choice belong to the owning recipe.
This samples the declared corpus; it does not establish an unbiased census of
prokaryotes, biological independence, or genome-payload identity.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.SelectNcrnaBenchmarkInputs "$@"
