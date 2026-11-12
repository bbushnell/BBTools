#!/bin/bash
# Report coverage from existing native sweep and reference metadata.
usage(){
cat <<'HELP'
Usage: dev/summarizencrnasweep.sh in=sweep.inputs.tsv changes=sweep.changes.tsv summary=sweep.summary.tsv reference=inputs.tsv loci=inventory.loci.tsv.gz members=members.tsv family=cobalamin policies=cap96,ceiling508 out=report

Joins saved paired losses/gains to original-input counts and reference support.
Reads metadata only; does not open FASTA/GFF payloads, run a caller or regrade.
All requested policies and baseline must cover the exact reference manifest.

baseline=baseline  Native policy used for paired lost/gained counts.
highlight=         Optional comma-separated reference-member names to flag.
out=prefix         Writes .loci.tsv, .inputs.tsv, .summary.tsv, .members.tsv.
taxonomy=          Optional saved full-corpus category table with input_id, fasta,
                   category. Adds .strata.tsv: first-gain/last-loss input counts
                   for every policy and each of the four saved categories.
                   Joins every on-panel ID/FASTA exactly once; never recomputes
                   rarity on the smaller panel. UNKNOWN and UNCLASSIFIED stay
                   separate and still count in aggregate detection coverage.
compare=left,right Optional exact scored-reference set comparison of two named
                   policies (including baseline), requiring taxonomy. Adds
                   .comparison.loci.tsv with original input/reference keys and
                   .comparison.summary.tsv with category/all losses, gains,
                   exclusive keys and subset flags. A net TP difference does
                   not imply set inclusion. This does not compare all GFF calls.
ow=f               Existing reports are protected by default; append unsupported.

Reports distinguish changed loci from distinct original inputs losing their
last detection or gaining their first. Input IDs are not collapsed into guessed
organisms. Member support is overlapping and nonadditive; marginal support stays
separate. Detection coverage is not sequence redundancy or subtype precision.
No scientific defaults are selected. One processing thread; -Xmx1g is a starting
heap choice for small sweep tables. File compression follows native suffixes.
HELP
}
if [[ $# == 0 || "$1" == -h || "$1" == --help ]]; then
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.SummarizeNcrnaSweep "$@"
