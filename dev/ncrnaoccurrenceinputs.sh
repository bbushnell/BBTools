#!/bin/bash
# Prepare native extraction inputs from a complete admitted occurrence table.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 7, 2026

Usage: dev/ncrnaoccurrenceinputs.sh in=occurrences.tsv roster=genomes.tsv out=new_directory expected=43013 flank=50

Prepares source-specific CutGff inputs, exact core coordinates and a distinct
learning ledger for admitted occurrences. No sequence extraction or filtering.
All occurrences must join the exact source genome, assembly and TaxID. Sources
with no admitted occurrences remain in the roster with zero loci.

in=        Complete NcrnaOccurrences occurrences.tsv.
roster=    Full eleven-column genome roster from NcrnaAnnotatedCorpus.
out=       Fresh directory beneath an existing parent; no overwrite/append.
expected=  Required exact occurrence count; this is a completeness check.
flank=50   Requested genomic flank per side, clipped to actual contig bounds.

Outputs: sources.config, gff/<source_index>.gff, core.gff, genomes.tsv,
occurrence_loci.tsv, coverage_split.tsv and summary.tsv.
Use sources.config with cutgff.sh to extract flanks with gccontig=caller, then
core.gff to derive exact cores from that FASTA. Preserve full FASTA descriptions.
The learning ledger separates FASTA record IDs from admitted occurrence IDs;
it does not invent historical annotation, partial or pseudogene labels.
One metadata worker; expected_context_bases is geometry, not measured extraction.

See docs/guides/NcrnaOccurrenceInputsGuide.txt for the complete contract.
Standard BBTools config=, outstream= and JVM options are supported.
Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter problems.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
while [ -h "$SCRIPT" ]; do
	DIR="$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$SCRIPT")"
	SCRIPT="$RESOLVED_SCRIPT"
	if [ "${SCRIPT#/}" = "$SCRIPT" ]; then
		RESOLVED_SCRIPT="$DIR/$SCRIPT"
		SCRIPT="$RESOLVED_SCRIPT"
	fi
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
if [ -f "$DIR/bbtools.jar" ]; then CP="$DIR/bbtools.jar"; else CP="$DIR/current/"; fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaOccurrenceInputs "$@"
