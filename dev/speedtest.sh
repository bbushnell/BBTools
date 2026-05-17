#!/bin/bash

usage(){
echo "
Written by Nahida
Last modified May 2026

Description:  ALU-constrained throughput benchmark for cardinality estimators.
Each thread processes one estimator at a time: create, add maxCard elements,
get cardinality, discard.  This keeps each estimator's working set in L1/L2
cache, measuring pure compute-bound insertion throughput.

Reports million adds/second for each estimator type.

Usage:  speedtest.sh buckets=2048 estimators=16384 card=40000000 t=8

Parameters:
buckets=2048      Bucket count for all estimator types.
estimators=16384  Total number of estimator instances to process.
card=40000000     Elements to add per estimator instance.
t=N               Threads (default: all available cores).
type=             Run only this type (default: run all types).

Types tested (in order):
  avll, exa, udll6, ull, dll4, ll6, hll4, hlll, htc4, htb

Notes:
  - This measures ALU-constrained throughput (1 estimator active per thread).
  - For bandwidth-constrained (cache pressure) throughput, use speedtest2.sh.
  - card should be large enough that the test runs >=10 seconds per type.
  - Output goes to stdout (TSV) and stderr (human-readable).

Java Parameters:
-Xmx            Heap size, e.g. -Xmx8g.

Contact: bbushnell@lbl.gov
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
	DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=2g" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.SpeedTest $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
