#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 22, 2025

Description:  Parses sequence files.
Reports bases and records.

Usage:  fastqscan.sh <file>

Input may be fastq, fasta, or sam, compressed or uncompressed.

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

EA="-da"
SIMD="--add-modules jdk.incubator.vector"
XMX="-Xmx256m"
XMS="-Xms256m"

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=1g" "--xms=256m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP stream.FastqScan $@"
	#echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
#setEnv "$@"
launch "$@"
