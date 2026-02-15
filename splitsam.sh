#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 9, 2015

Description:  Splits a sam file into three files:
Plus-mapped reads, Minus-mapped reads, and Unmapped.
If 'header' is the 5th argument, header lines will be included.

Usage:  splitsam <input> <plus output> <minus output> <unmapped output>

Input may be stdin or a sam file, raw or gzipped.
Outputs must be sam files, and may be gzipped.
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

	parseJavaArgs "--xmx=128m" "--xms=128m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.SplitSamFile $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
