#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 3, 2019

Description:  Shrinks or expands homopolymers.

Usage:  adjusthomopolymers.sh in=<input file> out=<output file> rate=<float>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
out=<file>      Primary output, or read 1 output.
out2=<file>     Read 2 output if reads are in two files.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
rate=0.0        0.1 will expand by 10%; -0.1 will shrink by 10%.

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

	parseJavaArgs "--xmx=1g" "--xms=1g" "--percent=30" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.AdjustHomopolymers $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
