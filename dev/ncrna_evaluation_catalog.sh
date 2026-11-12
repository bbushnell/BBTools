#!/bin/bash
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Raiden
Last modified November 8, 2026

Usage: dev/ncrna_evaluation_catalog.sh in=sources.tsv out=new_directory

Joins metadata from explicitly listed evaluation sources. Merges only identical
resolved filesystem paths; reports possible assembly aliases without merging them.
It does not read genome sequences, verify supplied hashes, or construct RNA truth.

in=              TSV: source, format, compartment, path, expected_records.
                 Optional sixth column digest_kind: file_bytes or unknown
                 (default). file_bytes declares hashes of the actual file bytes,
                 including compression, not normalized sequence hashes. Declare
                 it only from verified producer provenance.
out=             New directory; writes genomes.tsv, memberships.tsv, aliases.tsv,
                 sources.tsv and summary.tsv. Existing directories are never overwritten.
requirefiles=t   Fail after writing diagnostic outputs if any FASTA path is missing,
                 unreadable, empty or not a regular file. f permits an inventory
                 with unavailable paths; it does not certify those inputs.
resolvealiases=f Write alias_resolution.tsv when true. Reports a candidate only
                 when exact versioned assembly claims and recorded FASTA hashes
                 agree and exactly one path is available. Missing or conflicting
                 metadata and multiple live paths remain explicit. Report only:
                 memberships and requirefiles behavior do not change. Recorded
                 hashes are not checked against genome contents. sources.tsv
                 includes digest_kind in this mode. Candidate paths require
                 retained actual-file verification or new verification before
                 being accepted as scan inputs; this report is not that evidence.
ow=f             Required; no append, sampling or read limits.
-Xmx512m         Default heap for metadata, not genome sequence storage.

Source formats:
header     Tab-separated roster with a named fasta column. Optional columns:
           genome, assembly, taxid, compartment, gff and one of sha80,
           fasta_sha80 or selected_fasta_sha80. Extra columns are ignored.
config     Read direct in= entries from a caller configuration, one path per line.
           Nested config= and comma lists are rejected; no caller options execute.
taxonomy7  Headerless domain, phylum, taxid, genus, species, fasta, genome.

Source paths resolve relative to the input specifications; relative FASTA paths
resolve against the supplied source file's parent directory. Membership indices
count input records, not file
line numbers. GFF paths and digests are reported metadata, not verified bindings.
PRESENT means readable nonempty regular file, not valid FASTA or sequence identity.
Explicit versioned GCA/GCF/NC accessions produce alias candidates. Different
TaxIDs, names or accession versions never trigger an automatic merge.
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
if [ -f "$DIR/bbtools.jar" ]; then
	CP="$DIR/bbtools.jar"
else
	CP="$DIR/current/"
fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaEvaluationCatalog "$@"
