#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Eru
Last modified April 6, 2026

Description:  Measures DLC per-tier absolute error as a function of tier occupancy.
Uses LogLog6 (no minZeros, no microIndex, no tier promotion) for clean measurement.
Tier occupancy = number of buckets with NLZ >= tier, an integer from 0 to B.
Multithreaded: distributes DDL instances across threads.

Usage:  dlctieraccuracy.sh buckets=2048 ddls=100000 tier=3 threads=128

Parameters:
buckets=2048    Buckets per LL6 instance.
ddls=100000     Number of LL6 instances (distributed across threads).
maxmult=512     Stop when trueCard reaches buckets*maxmult.
tier=3          Which DLC tier to measure (0-63).
seed=1          Master seed for instance seed generation.
threads=1       Number of parallel simulation threads.
points=500      Number of log-spaced cardinality checkpoints.
type=ll6        Estimator type (usually ll6 for clean tier measurement).

Output columns (stdout, TSV):
Occupancy       Integer tier occupancy (0..buckets).
Count           Number of samples at this occupancy.
AvgAbsErr       Mean |DLC_tier_est - trueCard| / trueCard.
AvgSignedErr    Mean (DLC_tier_est - trueCard) / trueCard.

Java Parameters:
-Xmx            Memory limit.  -Xmx20g = 20 GB.
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.DLCTierAccuracy $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
