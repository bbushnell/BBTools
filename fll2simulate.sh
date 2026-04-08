#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Chloe von Einzbern-Bushnell
Last modified April 8, 2026

Description:  Simulates a single FLL2 (FutureLogLog 2-bit) word to build
per-tier correction factor tables.  Runs iters independent trials, each
feeding random hashes into one 16-bit FLL2 bank (6 buckets sharing a 4-bit
exponent) until the bank's local exponent exceeds maxTier.  Accumulates
(tier, history) statistics across all trials, combines order-equivalent
states, and prints tier average cardinalities and per-state multipliers.

Word layout: [15:12]=localExp [11:0]=6 x 2-bit future bitmap
Per-bucket:  LSB='seen floor hit'  MSB='seen floor+1 hit'
Promotion fires when all 6 LSBs set; MSBs shift to LSBs, localExp++.
Max cascade = 2.  IOT mode: hashes with delta>1 are ignored.

Usage:  fll2simulate.sh iters=10000 threads=8 maxTier=15

Parameters:
iters=10000     Number of simulation trials (more = better statistics).
threads=8       Number of parallel simulation threads.
maxTier=15      End each trial when localExp exceeds this value (0-15).
                Lower values run faster but miss high-tier statistics.

Java Parameters:
-Xmx            Override Java memory autodetection (e.g. -Xmx4g).
-eoom           Exit on out-of-memory exception (requires Java 8u92+).
-da             Disable assertions (faster, but skips correctness checks).

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.FLL2Simulator $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
