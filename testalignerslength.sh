#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Neptune
Last modified March 23, 2026

Description:  Benchmark search space vs sequence length at fixed ANI.
Generates random sequence pairs at each length, runs Glocal as truth,
then each heuristic aligner. Multithreaded with loop tracking.

Usage:
testalignerslength.sh ani=75 samples=100 threads=64 lengths=64,128,256,512

Parameters:
ani             Target ANI in percent (default 75).
samples         Number of random pairs per length (default 100).
threads         Number of parallel threads (default all available).
lengths         Comma-separated list of sequence lengths.
subsonly        Substitutions only, no indels (default false).
equalrates      Equal S/D/I rates at 33/33/33 (default false).
seed            Random seed for reproducibility (default 54321).

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
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

	parseJavaArgs "--xmx=2000m" "--xms=2000m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP idaligner.TestAlignerLength $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
