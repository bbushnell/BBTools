#!/bin/bash

usage(){
echo "
Written by Chloe
Last modified October 18, 2025

Description:  Converts BAM (Binary Alignment/Map) files to SAM 
(Sequence Alignment/Map) text format. Reads BGZF-compressed BAM files 
and outputs tab-delimited SAM format.

Usage:  bamlinestreamer.sh <input.bam> <output.sam>

Standard parameters:
in=<file>        Input BAM file (first positional argument).
out=<file>       Output SAM file (second positional argument).

Java Parameters:
-Xmx             This will set Java's memory usage, overriding autodetection.
                 -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                 specify 200 megs. The max is typically 85% of physical memory.
-eoom            This flag will cause the process to exit if an out-of-memory
                 exception occurs.  Requires Java 8u92+.
-da              Disable assertions.

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

	parseJavaArgs "--xmx=8g" "--xms=8g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP stream.bam.Bam2Sam $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
