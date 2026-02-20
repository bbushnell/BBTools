#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 1, 2016

Description:  Generates a kmer uniqueness histogram, binned by file position.
There are 3 columns for single reads, 6 columns for paired:
count        number of reads or pairs processed
r1_first     percent unique 1st kmer of read 1
r1_rand      percent unique random kmer of read 1
r2_first     percent unique 1st kmer of read 2
r2_rand      percent unique random kmer of read 2
pair         percent unique concatenated kmer from read 1 and 2

Please read bbmap/docs/guides/CalcUniquenessGuide.txt for more information.

Usage:	bbcountunique.sh in=<input> out=<output>

Input parameters:
in2=null            Second input file for paired reads
interleaved=auto    Set true/false to override autodetection of the input file as paired interleaved.
samplerate=1        Set to below 1 to sample a fraction of input reads.
reads=-1            Only process this number of reads, then quit (-1 means all)

Output parameters:
out=<file>          File for output stats

Processing parameters:
k=25                Kmer length (range 1-31).
interval=25000      Print one line to the histogram per this many reads.
cumulative=f        Show cumulative numbers rather than per-interval numbers.
percent=t           Show percentages of unique reads.
count=f             Show raw counts of unique reads.
printlastbin=f      (plb) Print a line for the final undersized bin.
minprob=0           Ignore kmers with a probability of correctness below this (based on q-scores).

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.CalcUniqueness $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
