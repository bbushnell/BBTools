#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 14, 2025

Description:  Aligns a query sequence to a reference using ScrabbleAligner.
The sequences can be any characters, but N is a special case.
Outputs the identity, rstart, and rstop positions.
Optionally prints a state space exploration map.
This map can be fed to visualizealignment.sh to make an image.

Usage:
scrabblealigner.sh <query> <ref>
scrabblealigner.sh <query> <ref> <map>
scrabblealigner.sh <query> <ref> <map> <iterations>

Parameters:
query           A literal nucleotide sequence or fasta file.
ref             A literal nucleotide sequence or fasta file.
map             Optional output text file for matrix score space.
                Set to null for benchmarking with no visualization.
iterations      Optional integer for benchmarking multiple iterations.

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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.ScrabbleAligner $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"