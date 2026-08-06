#!/bin/bash

usage(){
echo "
Written by Eru
Last modified July 31, 2026

Description:  Blastp-style protein similarity search (Phase-1 MVP).
Reads a query protein FASTA and a database protein FASTA and writes
BLAST-tab outfmt 6 TSV results.  Seeds candidates with exact (or amino8
reduced-alphabet) k-mers, then scores with BLOSUM62 affine-gap local
alignment (gap-open 11, gap-extend 1).

Usage:  proteinsearch.sh query=<query.faa> db=<database.faa> out=<hits.tsv>

Required parameters:
query=<file>    Query protein FASTA (amino acids).
db=<file>       Database protein FASTA (amino acids).

Optional parameters (and their defaults):
out=stdout      Output TSV (outfmt 6).  A .meta sidecar is written beside it.
evalue=10       E-value significance cutoff.
minid=0         Minimum percent identity to report.
minscore=0      Minimum raw BLOSUM62 score to report.
k=5             Seed k-mer length.
reduced=f       Use amino8 reduced-alphabet seeds (more sensitive).
mts=            max-target-seqs: cap distinct targets per query.
ow=t            (overwrite) Overwrite existing output.

Output columns (tab-separated):
query target pident length mismatch gapopen qstart qend tstart tend evalue bitscore

Notes:
Coordinates are 1-based inclusive.  The bitscore is rigorous (gapped
BLOSUM62 11/1); the E-value is approximate (Karlin-Altschul edge-length
correction omitted) and is flagged in the .meta sidecar.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
-eoom           Exit if an out-of-memory exception occurs.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"

}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=2g" "--xms=2g" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.ProteinSearch $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
