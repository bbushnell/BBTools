#!/bin/bash

usage(){
echo "
Shell script written by Brian Bushnell
Java code written by Claude.
Last modified May 4, 2025

Description:  Converts a text exploration map from some aligners to an image.
Supports Quantum, Banded, Drifting, Glocal, WaveFront, and MSA9.

Usage:
visualizealignment.sh <map>
or
visualizealignment.sh <map> <image>

Parameters:
map             Text file of score-space from an aligner.
image           Output name, context sensitive; supports png, bmp, jpg.
                Image name is optional; if absent, .txt will be replaced
                by .png in the input filename.

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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.VisualizationConverter $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
