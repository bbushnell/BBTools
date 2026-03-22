#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Neptune
Last modified March 22, 2026

Description:  Batch benchmark for aligners using random sequences.
Pre-generates all sequence pairs, runs Glocal first to establish truth,
then runs each other aligner. Uses multithreaded work-stealing.

Usage:
testalignersbatch.sh length=40000 samples=100 threads=64 subsonly=t ani=100,99,95,90

Parameters:
length          Sequence length in bp (default 40000).
samples         Number of random pairs per ANI level (default 10).
threads         Number of parallel threads (default all available).
subsonly        Substitutions only, no indels (default false).
equalrates      Equal S/D/I rates at 33/33/33 (default false).
ani             Comma-separated list of design ANI values.
seed            Random seed for reproducibility (default 12345).

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP idaligner.TestAlignerBatch $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
