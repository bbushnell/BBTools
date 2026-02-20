#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 12, 2025

Description:  Fast lightweight scanner that parses newlines.
Supports raw, gzip, bgzip, and bz2 compression, and any text filetype.

Usage:  filescan.sh <file> <threads>
e.g.
filescan.sh contigs.fasta
filescan.sh reads.fq.gz
filescan.sh reads.fq 2

Bgzipped input processing is multithreaded and much faster than regular gzip.
SIMD support is autodetected and can be disabled with the flag simd=f.

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

EA="-ea"
SIMD="--add-modules jdk.incubator.vector"
XMX="-Xmx256m"
XMS="-Xms128m"

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=256m" "--xms=128m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP stream.FileScanMT $@"
	#echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
#setEnv "$@"
launch "$@"
