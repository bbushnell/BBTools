#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 20, 2015

Description:  Prints the most common kmers in each sequence.
This is intended for short kmers only!

Usage:  commonkmers.sh in=<file> out=<file>

Parameters:
k=2             Kmer length, 0-12.
display=3       Print this many kmers per sequence.
count=f         Print the kmer counts as well.

ow=f            (overwrite) Overwrites files that already exist.
app=f           (append) Append to files that already exist.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).
qin=auto        ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
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

	parseJavaArgs "--xmx=800m" "--xms=800m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.SmallKmerFrequency $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
