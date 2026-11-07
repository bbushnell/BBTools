#!/bin/bash
# Check complete extracted RNA context/core populations against their ledger.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 7, 2026

Usage: ncrnaoccurrencecheck.sh in=occurrence_loci.tsv contexts=contexts.fa cores=cores.fa expected=43013 out=summary.tsv

Checks exact record-ID sets, source tags, actual context/core lengths, declared
sense flanks and byte-for-byte agreement of every core with its actual context.
Context GC must be finite in [0,1] and identical for records from one source
contig. GC comes from the native caller-compatible extractor; this does not
independently recount whole genomes or establish biological annotation truth.

in=        Fifteen-column ledger from ncrnaoccurrenceinputs.sh.
contexts=  Actual flanked FASTA with source/lflank/rflank/contig_gc headers.
cores=     Actual exact-core FASTA with the same record IDs and source tags.
expected=  Required full population count. Missing/duplicate/foreign IDs fail.
out=       Fresh summary TSV, opened only after all checks pass.

Both short-RNA FASTAs are retained through native sequence readers before
checking. Memory scales with total sequence/header size; this is not intended
for loading whole genomes. No sampling, exclusions, sequence normalization or
biological thresholds. Preserve full headers. Standard BBTools JVM options apply.
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
DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
if [ -f "$DIR/bbtools.jar" ]; then CP="$DIR/bbtools.jar"; else CP="$DIR/current/"; fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaOccurrenceCheck "$@"
