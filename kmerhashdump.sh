#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 5, 2026

Description:  Reads a fasta/fastq file, encodes every kmer as a 2-bit-packed
long using the standard shift-mask approach, and anonymizes it with
Tools.hash64shift.  Dumps the resulting hashcodes one per line (A48-encoded),
so downstream tools can see the statistical structure of real data (repeats,
duplicates, complexity) without ever handling or reasoning about the
underlying sequence.  Built for testing cardinality estimators against
realistic hash streams.

Usage:  kmerhashdump.sh in=<file> out=<file>

Input may be fasta or fastq, compressed or uncompressed.
Output is a text file, one A48-encoded hash per line.

Parameters:
in=<file>       Input sequence file.
out=<file>      Output file of hashcodes, one per line.
k=31            Kmer length (1-31; must fit in a 2-bit-packed long).
reads=-1        Only process this number of reads, then quit (-1 means all).
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
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

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.KmerHashDump $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
