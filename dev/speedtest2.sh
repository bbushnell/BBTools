#!/bin/bash

usage(){
echo "
Written by Nahida
Last modified May 2026

Description:  Bandwidth-constrained throughput benchmark for cardinality estimators.
Each thread maintains N simultaneous estimators and cycles through all of them on
every add, creating cache pressure that simulates real-world deployment with many
active sketches (e.g., per-reference or per-kmer-group tracking).

Reports million adds/second for each estimator type.

Usage:  speedtest2.sh mem=2k sim=128 card=40000000 t=8

Parameters:
mem=2k          Memory budget per estimator (e.g. 2k, 1024).  Each type
                gets the maximum bucket count that fits in this budget.
                Overrides 'buckets' when set.
buckets=2048    Bucket count for all estimator types (ignored when mem is set).
sim=128         Simultaneous estimators per thread.
card=40000000   Elements to add per estimator instance.
t=N             Threads (default: all available cores).
type=           Run only this type (default: run all types).

Types tested (in order):
  avll, exa, udll6, ull, dll4, ll6, hll4, hlll, htc4, htb

Notes:
  - This measures bandwidth-constrained throughput (N estimators active per thread).
  - For ALU-constrained (single estimator) throughput, use speedtest.sh.
  - sim controls cache pressure: sim=1 is ~speedtest.sh, sim=4096 is extreme.
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.SpeedTest2 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
