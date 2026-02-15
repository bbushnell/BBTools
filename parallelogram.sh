#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 4, 2025

Description:  Converts a parallelogram-shaped alignment visualization to a rectangle.
This tool transforms the output from CrossCutAligner so it can be properly
visualized by visualizealignment.sh. The transformation shifts coordinates
to create a rectangular matrix from the parallelogram pattern.

Usage:
parallelogram.sh <input_map> <output_map>

Parameters:
input_map       Input text file containing parallelogram-shaped matrix data.
output_map      Output text file with rectangular matrix data.

Example workflow:
crosscutaligner.sh ATCGATCG GCATGCTA map1.txt
parallelogram.sh map1.txt map2.txt
visualizealignment.sh map2.txt alignment.png

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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.Parallelogram $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"