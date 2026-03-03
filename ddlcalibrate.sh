#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Chloe
Last modified March 2, 2026

Description:  Calibrates DynamicDemiLog cardinality estimators by feeding random
longs to an array of DDL instances with varied seeds, tracking true cardinality
via ground truth, and reporting per-estimator accuracy statistics at logarithmically
spaced cardinality checkpoints.  Multithreaded: N threads each run an independent
simulation, results are merged live for File 1 output.  At end, a second file
(out2=) captures the occupancy histogram for per-occupancy compensation curve fitting.

Usage:  ddlcalibrate.sh ddls=1000 buckets=2048

Parameters:
ddls=1000       Number of DDL instances with varied seeds (distributed across threads).
buckets=2048    Buckets per DDL instance. Must be a multiple of 256.
maxmult=10      Stop when trueCard reaches buckets*maxmult.
reportfrac=0.01 Report every this fraction of current cardinality (cascading:
                gives one row per add at low cardinality, one per 10 at 1000, etc).
seed=12345      Master seed for DDL seed generation.
valseed=42      Master seed for per-thread value generation.
threads=1       Number of parallel simulation threads.
out2=           Output file for occupancy histogram (File 2).  Not written if omitted.

File 1 output columns (stdout, live):
TrueCard        Ground-truth distinct count.
Occupancy       Fraction of DDL buckets with any data (averaged across all DDLs).
*_err           Signed relative error: (estimate - true) / true.  Zero = unbiased.
*_abs           Mean absolute relative error.  Zero = perfectly accurate.
*_std           Population stdev of relative error across DDL instances.  Zero = no variance.

File 2 output columns (out2=, end of run):
Slot            Occupancy slot index 0..256 (256 = 100% full).
Occupancy       Exact occupancy at this slot (slot/256).
AvgTrueCard     Average true cardinality across all samples at this slot.
Samples         Number of DDL snapshots recorded at this slot.
*_raw           Average raw (uncompensated) estimate at this slot.
*_cf            Compensation factor: AvgTrueCard / *_raw.  Apply at runtime for correction.

Estimators reported (11 total):
Mean            Arithmetic mean proxy, with bucket correction factor.
HMean           HLL-style harmonic mean, with correction.
GMean           Geometric mean proxy, with correction.
RMean           Root-mean-square proxy, with correction.
MWA             Median-weighted average proxy (odd-list bug pending fix).
MedianCorr      Median proxy, with correction.
MedianLeg       Median proxy, legacy (no correction factor).
EstSum          Log-sum estimator (div * exp(mean of log estimates)).
LCHybrid        LinearCounting falling back to Mean at full occupancy.
LCTrue          Pure LinearCounting (accurate only at low occupancy).
Blended         Sigmoid blend of LCTrue and Mean; current cardinality() output.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.DDLCalibrationDriver $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
