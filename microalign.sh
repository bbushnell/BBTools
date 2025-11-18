#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 21, 2025

Description:  Wrapper for MicroAligner.
Can align reads to a small, single-contig reference like PhiX.
Probably faster than BBMap.  Produces most of the same histograms,
like idhist, mhist, etc.
Not currently designed for reference with multiple sequences,
or duplicate kmers of length used for indexing.

Usage:  microalign.sh in=<input file> out=<output file> ref=<reference>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
out=<file>      Primary output, or read 1 output.
out2=<file>     Read 2 output if reads are in two files.
outu=<file>     Optional unmapped read output.
outu2=<file>    Optional unmapped read 2 output.

Processing parameters:
k=17            Main kmer length.
k2=13           Sub-kmer length for paired reads only.
minid=0.66      Minimum alignment identity.
minid2=0.56     Minimum alignment identity if the mate is mapped.
mm=1            Middle mask length; the index uses gapped kmers.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
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
	SCRIPT="$0"
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.MicroWrapper $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
