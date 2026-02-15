#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 18, 2025

Description:  Aligns a query sequence to a reference using multiple aligners.
Outputs the identity, rstart and rstop positions, time, and #loops.

Usage:
testaligners.sh <query> <ref>
testaligners.sh <query> <ref> <iterations> <threads> <simd>

Parameters:
query           A literal nucleotide sequence or fasta file.
ref             A literal nucleotide sequence or fasta file.
iterations      Optional integer for benchmarking multiple iterations.
threads         Number of parallel instances to use.
simd            Enable SIMD operations; requires AVX-256 and Java 17+.

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

	parseJavaArgs "--xmx=2000m" "--xms=2000m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP idaligner.Test $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
