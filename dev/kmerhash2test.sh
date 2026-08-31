#!/bin/bash

usage(){
echo "
Written by Noelle
Last modified August 2026

Description:  Regression and distribution checks for Tools.splitMix64,
ukmer.Kmer.xor2(), and resizable/fixed HashArrayH1D storage, including
reverse-complement invariance, lazy invalidation, restoration, collisions,
bit balance, resize conservation, ownership, deletion, and loud capacity failure.

Usage:  kmerhash2test.sh

Contact: bbushnell@lbl.gov
"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
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
	DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"
	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch(){
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ukmer.KmerHash2Test"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
