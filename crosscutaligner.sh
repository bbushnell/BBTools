#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 4, 2025

Description:  Aligns a query sequence to a reference using CrossCutAligner.
This fully explores the matrix using 4 arrays of roughly length reflen.
The sequences can be any characters, but N is a special case.
Outputs the identity, rstart, and rstop positions.
CrossCut is a nontraditional aligner that fills antidiagonals,
incurring zero data dependencies between loops.  This allows
perfect SIMD vectorization.

Usage:
crosscutaligner.sh <query> <ref>
crosscutaligner.sh <query> <ref> <map>
crosscutaligner.sh <query> <ref> <map> <iterations> <simd>

Parameters:
query           A literal nucleotide sequence or fasta file.
ref             A literal nucleotide sequence or fasta file.
map             Optional output text file for matrix score space.
                Set to null for benchmarking with no visualization.
                This has not yet been tested and will produce unknown results.
iterations      Optional integer for benchmarking multiple iterations.
simd            Use vector instructions.

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP idaligner.CrossCutAligner $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
