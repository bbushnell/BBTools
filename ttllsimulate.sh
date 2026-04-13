#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Ady
Last modified April 2026

Description:  Simulates a single TTLL (TwinTailLogLog) word to build a
per-tier state table and per-state correction factors.  Each word is an
8-bit register: [7:4]=4-bit exponent, [3:2]=history1, [1:0]=history0.
Each trial feeds random hashes into one bucket until its local exponent
exceeds maxTier.  The simulator accumulates (tier, combined_history)
statistics across trials, smooths sparse states, and emits tier averages
plus per-state multipliers for CF table construction.

Usage:  ttllsimulate.sh iters=10000 threads=8 maxTier=14

Parameters:
iters=10000     Number of simulation trials (more = better statistics).
threads=8       Number of parallel simulation threads.
maxTier=14      End each trial when localExp exceeds this value (0-14).
minObs=100      Merge states with fewer than this many observations.
out=            Optional output file for state table (default: stdout).
table=          Optional input state table for CV reference.
avg=lin         Averaging mode: lin | geo | harm | blend.

Java Parameters:
-Xmx            Override Java memory autodetection (e.g. -Xmx4g).
-eoom           Exit on out-of-memory exception (requires Java 8u92+).
-da             Disable assertions.

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.TTLLSimulator $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
