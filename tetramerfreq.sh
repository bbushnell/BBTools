#!/bin/bash

usage(){
echo "
Written by Shijie Yao and Brian Bushnell
Last modified April 25, 2025

Description: DNA Tetramer analysis.
DNA tetramers are counted for each sub-sequence of window size in the sequence.
The window slides along the sequence by the step length.
Sub-sequence shorter than the window size is ignored. Tetramers containing N are ignored.

Usage: tetramerfreq.sh in=<input file> out=<output file> step=500 window=2000

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       DNA sequence input file
out=<file>      Output file name
step/s=INT      Step size (default 500)
window/w=INT    Window size (default 2kb); <=0 turns windowing off (e.g. short reads)
short=T/F       Print lines for sequences shorter than window (default F)
k=INT           Kmer length (default 4)
gc              Print a GC column in the output.
float           Output kmer frequencies instead of counts.
comp            Output GC-compensated kmer frequencies.

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

	parseJavaArgs "--xmx=1000m" "--xms=1000m" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.TetramerFrequencies $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
