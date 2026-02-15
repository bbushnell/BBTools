#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 5, 2020

Description:  Aligns a query sequence to reference sequences.
Outputs the best matching position per reference sequence.
If there are multiple queries, only the best-matching query will be used.
MSA in this context stands for MultiStateAligner, not Multiple Sequence Alignment.

Usage:
msa.sh in=<file> out=<file> literal=<literal,literal,...>
or
msa.sh in=<file> out=<file> ref=<lfile>

Parameters:
in=<file>       File containing reads.
out=<file>      Sam output file.
literal=        A sequence of bases to match, or a comma-delimited list.
ref=<file>      A fasta file of bases to match.  Please set either ref
                or literal, not both.
rcomp=t         Also look for reverse-complements of the sequences.
addr=f          Add r_ prefix to reverse-complemented alignments.
replicate=f     Make copies of sequences with undefined bases for every
                possible combination.  For example, ATN would expand to
                ATA, ATC, ATG, and ATT.
cutoff=0        Ignore alignments with identity below this (range 0-1).
swap=f          Swap the reference and query; e.g., report read alignments
                to the reference instead of reference alignments to the reads.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding automatic
                memory detection. -Xmx20g will specify
                20 gigs of RAM, and -Xmx200m will specify 200 megs.
                The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
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
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=2g" "--xms=2g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.FindPrimers $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
