#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Chloe
Last modified March 16, 2026

Description:  Tests cardinality estimator accuracy on low-complexity datasets
with bounded cardinality and repeated values.  Draws with replacement from a
fixed array of unique values, biased toward lower indices via min(rand, rand)
to simulate skewed frequency distributions.  Estimates are recorded after
every add to capture behavior while 'parked' at a given cardinality.

Usage:  lowcomplexcalibrate.sh card=5000 ddls=128 type=dll3

Parameters:
card=5000       Maximum true cardinality (number of unique values).
ddls=128        Number of estimator instances with varied seeds.
buckets=2048    Buckets per estimator. Must be a multiple of 256.
iter=0          Iterations as a multiplier of cardinality. 0=stop on saturation.
type=dll4       Estimator type: ddl, ddl2, ddl8, dll2, dll3, dll3v2, dll4.
threads=4       Number of parallel threads.
seed=1          Master seed for value array and estimator seeds.
cf=t            Enable/disable correction factors.
cardcf=t        Enable/disable cardinality-based correction factors.

Java Parameters:
-Xmx            Set Java's memory usage (e.g. -Xmx4g).
-eoom           Exit on out-of-memory exception. Requires Java 8u92+.
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.LowComplexityCalibrationDriver $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
