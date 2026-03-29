#!/bin/bash

usage(){
echo "
Written by Eru
Last modified March 29, 2026

Description:  Trains an LC history correction table for history-aware Linear Counting.
Simulates HyperLogLog with N-bit history, recording per-bucket state distributions
at each occupancy level.  Each thread runs independent trials with its own RNG and
accumulators; results are merged after all threads finish.

Usage:  trainLCHist.sh buckets=2048 t=4 hbits=2 iters=1m out=table.tsv.gz

Parameters:
buckets=2048    Number of buckets (must be power of 2).
hbits=2         Number of history bits (1, 2, or 3).
iters=100k      Number of simulation trials (supports K/M/G suffixes).
t=1             Number of parallel threads.
seed=42         Master seed for deterministic results.
out=stdout.txt  Output file.  Use .gz extension for gzip compression.

Java Parameters:
-Xmx            Memory limit.  -Xmx1g is usually sufficient.
-eoom           Exit on out-of-memory.
-da             Disable assertions.
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
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=1g" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.LCHistTrainer $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
